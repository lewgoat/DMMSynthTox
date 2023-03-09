#import tensorflow as tf
from pandas import read_csv, DataFrame as df
#from keras import *
#from keras.layers import *
#from mol2vec import features
from utils import featurize
from gensim.models import Word2Vec
import shutil

TRAINING_SET_SIZE = 1000
MOL_TO_VEC_DIM_PREFIX = "mol2vec-"

def smiles_to_vec():
    featurize('datatable.smi', 'vectorised.csv', 'model_300dim.pkl', r=300)
    
    vectors = read_csv("vectorised.csv", usecols=lambda x: x.startswith(MOL_TO_VEC_DIM_PREFIX))
    return vectors

def main():

    smiles_to_vec()

    """
    net = Sequential()
    net.add(Input(shape=(16,TRAINING_SET_SIZE)))
    net.add(Dense(128, activation="relu", name="hidden"))
    net.add(Dense(1, name="output"))

    net.predict()
    """

if __name__ == "__main__":
    main()