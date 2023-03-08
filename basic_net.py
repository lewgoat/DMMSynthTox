#import tensorflow as tf
from pandas import read_csv, DataFrame as df
#from keras import *
#from keras.layers import *
from mol2vec import features
from gensim.models import Word2Vec
import shutil

TRAINING_SET_SIZE = 1000

def smiles_to_vec():

    data = read_csv("datatable.csv", usecols=["Smiles", "ID"], index_col=False).iloc[:, ::-1]
    print(data.head(10))
    data_tab = open("datatable.smi", "a")
    data.to_csv(data_tab, sep='\t', lineterminator='\n', index=False)
    data_tab.close()

    features.featurize('datatable.smi', 'vectorised.csv', 'model_300dim.pkl', r=300)

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