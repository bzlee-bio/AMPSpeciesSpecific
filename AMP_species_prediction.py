from Bio import SeqIO
import pandas as pd
import tensorflow as tf
import numpy as np
import pickle
import os
import argparse
import copy
import re
import math

### Data loader
class w2v:
    def __init__(self):
        self.vocab = {' ': 0,
                     'A': 1,
                     'C': 2,
                     'D': 3,
                     'E': 4,
                     'F': 5,
                     'G': 6,
                     'H': 7,
                     'I': 8,
                     'K': 9,
                     'L': 10,
                     'M': 11,
                     'N': 12,
                     'P': 13,
                     'Q': 14,
                     'R': 15,
                     'S': 16,
                     'T': 17,
                     'V': 18,
                     'W': 19,
                     'Y': 20}
    def w2vec(self, data):
        vec_data = np.zeros((len(data),100))
        lab = []
        rm_list = []
        seqs = []
        
        for i, l in enumerate(data):
            l[1] = l[1].replace(' ','').replace('*','')
            try:
                vec_data[i, :len(l[1])]=np.array([self.vocab[AA.upper()] for AA in list(l[1])])
                lab.append(l[0])
                seqs.append(l[1])
            except:
                print('Except!')
                print(i, l[0], l[1])
                rm_list.append(i)
        for r_i in rm_list:
            vec_data = np.delete(vec_data,r_i, 0)
        return vec_data, lab, seqs



        
### Model generator
class Model(tf.keras.Model):
    def __init__(self, params):
        super(Model, self).__init__()
        self._embedding = tf.keras.layers.Embedding(input_dim = params['EMB_input_dim'],output_dim=params['EMB_output_dim'],\
                                        mask_zero=True, name='EMBED')
        ## Feature extracting layer
        self._FE_LSTM = [tf.keras.layers.Bidirectional(tf.keras.layers.LSTM(units, return_sequences=True,\
                                        dropout=params['dropout'],\
                                        kernel_regularizer=tf.keras.regularizers.l1(0.01),\
                                        activity_regularizer=tf.keras.regularizers.l2(0.01)), name='FE_LSTM_'+str(i)) \
                            for i, units in enumerate(params['FE_LSTM_unit'])]

        ## Specific targeting layer
        self._ST_Layers = {}
        for tgt in params['tgt_list']:
            self._ST_Layers[tgt] = [tf.keras.layers.Bidirectional(tf.keras.layers.LSTM(units, return_sequences=True,\
                                        dropout=params['dropout'],\
                                        kernel_regularizer=tf.keras.regularizers.l1(0.01),\
                                        activity_regularizer=tf.keras.regularizers.l2(0.01)), name='FE_LSTM_'+str(i)+'_'+ tgt.replace("/",'').replace(' ','')) \
                                    for i, units in enumerate(params['ST_LSTM_unit'])]
            self._ST_Layers[tgt]+=[tf.keras.layers.Flatten()]
            
            for i, units in enumerate(params['DENSE_unit']):
                if i == (params['DENSE_unit'].__len__()-1):
                    self._ST_Layers[tgt]+=[tf.keras.layers.Dense(units, name='DENSE_'+str(i)+'_'+ tgt.replace("/",'').replace(' ',''))]
                else:
                    self._ST_Layers[tgt]+=[tf.keras.layers.Dense(units, name='DENSE_'+str(i)+'_'+ tgt.replace("/",'').replace(' ',''), \
                                                   activation='relu')]
        self.params = params
        
    def model_build(self):
        inp = tf.keras.Input(shape=(100))
        x = self._embedding(inp)
        for FE_LSTM in self._FE_LSTM:
            x = FE_LSTM(x)
        FE_out = x
        tgt_output = {}
        tgt_model = {}
        for tgt in self.params['tgt_list']:
            for i, ST_Layers in enumerate(self._ST_Layers[tgt]):
                if i==0:
                    y = ST_Layers(FE_out)
                else:
                    y = ST_Layers(y)
            tgt_output[tgt] = y
    
            tgt_model[tgt] = tf.keras.Model(inp, tgt_output[tgt])
#         print(list(tgt_model.values()))
        tot_output = tf.keras.layers.Concatenate(name='Total_output')(list(tgt_output.values()))
        tot_model = tf.keras.Model(inp, tot_output)
        return tgt_model, tot_model
                
                
            
### Model trainer
class model_train:
    def __init__(self, params):
        M = Model(params)
        self.tgt_model, self.tot_model = M.model_build()
        print(self.tot_model.summary())
        self.loss_fn = tf.keras.losses.BinaryCrossentropy(from_logits=True)
        self.optimizer = tf.keras.optimizers.Adam(params['Learning_rate'])
        self.params = params
        
  
        
    def prediction_prob_out(self, data):
        return tf.keras.activations.sigmoid(self.tot_model.predict(data))

    def prediction_raw_out(self, data):
        return self.tot_model.predict(data)
            
    def model_loader(self, load_file):
        load_status = self.tot_model.load_weights(load_file)
        print(load_status.assert_consumed())
            
            

parser = argparse.ArgumentParser(description="Prediction of antimicrobial activity for B. subtilis, E. coli, P. aeruginosa, S. aureus, S. epidermidis")
parser.add_argument('--windowsize', help='Size of window for sliding window. Without input of window size, the prediction performed for sequences length <= 50',default=None)
parser.add_argument('--step', help='Step size of sliding window (Default: 1)',default=1)
parser.add_argument('--fasta', help='Input fasta file')
parser.add_argument('--out', help='Output file')    

args = parser.parse_args()
fastaFile = args.fasta

if args.windowsize is not None:
    windowSize = int(args.windowsize)
else:
    windowSize = args.windowsize
if args.step is not None:
    step = int(args.step)
else:
    step = args.step
outFile = args.out
print(f'Sequence load..')
loadSeqs = []
if windowSize is None:
    for seq in SeqIO.parse(fastaFile, 'fasta'):
        if len(str(seq.seq))<=50:
            loadSeqs.append([seq.id, str(seq.seq)])
else:
    for seq in SeqIO.parse(fastaFile, 'fasta'):
        if len(str(seq.seq))<=50:
            loadSeqs.append([seq.id, str(seq.seq)])
        else:
            for i in range(math.ceil((len(str(seq.seq))-windowSize)/step)+1):
                
                startIdx = int(i*step)
                stopIdx = int(i*step+windowSize)
                loadSeqs.append([f'{seq.id}_window_{i}', str(seq.seq)[startIdx:stopIdx]])
            
            



w2vec = w2v()
embSeqs, ids, seqs = w2vec.w2vec(loadSeqs)
idsLen = len(ids)
print(f'Total {idsLen} sequences')
print(seqs)

modelParams = f'./modelWeights/params.pickle'
modelFile = f'./modelWeights/modelWeight'

with open(modelParams, 'rb') as f:
    params = pickle.load(f)

mt = model_train(params)
mt.model_loader(modelFile)
predRes = mt.prediction_prob_out(embSeqs).numpy()
predRes = pd.DataFrame(predRes, columns=['B. subtilis', 'E. coli','P. aeruginosa','S. aureus','S. epidermidis'])
idsSeries = pd.Series(ids)
seqSeries = pd.Series(seqs)
predRes.insert(0,'ID', idsSeries)
predRes.insert(1,'Sequence', seqSeries)                                 
predRes.to_csv(outFile, index=None)
