#import tensorflow as tf
from pandas import read_csv, DataFrame as df
from keras import Sequential, optimizers
from keras.layers import Dense
import numpy as np
import matplotlib.pyplot as plt
#from mol2vec import features
from utils import featurize
from gensim.models import Word2Vec
import shutil

NO_OF_FEATURES = 100
TRAINING_SET_SIZE = 1000
VALIDATION_PERCENT = 10
MOL_TO_VEC_DIM_PREFIX = "mol2vec-"

def smiles_to_vec():
    #featurize('datatable.smi', 'vectorised.csv', 'model_300dim.pkl', r=300)
    
    vectors = read_csv("vectorised.csv", usecols=lambda x: x.isnumeric())
    return vectors

def main():

    vec_frame = smiles_to_vec()
    inputs = vec_frame.to_numpy()
    data_size = np.shape(inputs)[0]
    valid_size = int((VALIDATION_PERCENT / 100) * data_size)
    train_size = data_size - valid_size
    ld_frame = read_csv("datatable.smi", usecols=["LD50"])
    ld = ld_frame.to_numpy()

    net = Sequential()
    net.add(Dense(128, input_dim=NO_OF_FEATURES, activation="relu", name="hidden1"))
    net.add(Dense(128, activation="relu", name="hidden2"))
    net.add(Dense(1, name="output"))

    net.compile(optimizer="rmsprop", loss='mean_squared_logarithmic_error')
    history = net.fit(
        inputs,
        ld,
        epochs=300,
        verbose=0,
        validation_split=0.1)

    plot_loss(history)

    #net.predict()
    

def plot_loss(history):
  plt.plot(history.history['loss'], label='loss')
  plt.plot(history.history['val_loss'], label='val_loss')
  plt.ylim([0, 2])
  plt.xlabel('Epoch')
  plt.ylabel('Error [MPG]')
  plt.legend()
  plt.grid(True)
  plt.show()

if __name__ == "__main__":
    main()