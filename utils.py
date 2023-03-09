import numpy as np
import pandas as pd
from gensim.models import word2vec
from rdkit import Chem
from rdkit.Chem import PandasTools
from mol2vec.features import MolSentence, mol2alt_sentence

def sentences2vec(sentences, model, unseen=None):
    """Generate vectors for each sentence (list) in a list of sentences. Vector is simply a
    sum of vectors for individual words.
    
    Parameters
    ----------
    sentences : list, array
        List with sentences
    model : word2vec.Word2Vec
        Gensim word2vec model
    unseen : None, str
        Keyword for unseen words. If None, those words are skipped.
        https://stats.stackexchange.com/questions/163005/how-to-set-the-dictionary-for-text-analysis-using-neural-networks/163032#163032

    Returns
    -------
    np.array
    """
    keys = set(model.wv.key_to_index)
    vec = []
    
    if unseen:
        unseen_vec = model.wv.get_vector(unseen)

    for sentence in sentences:
        if unseen:
            vec.append(sum([model.wv.get_vector(y) if y in set(sentence) & keys
                       else unseen_vec for y in sentence]))
        else:
            vec.append(sum([model.wv.get_vector(y) for y in sentence 
                            if y in set(sentence) & keys]))

    filtered_vec = [v if isinstance(v, np.ndarray) else np.zeros(100) for v in vec]
    return np.array(filtered_vec)

def featurize(in_file, out_file, model_path, r, uncommon=None):
    """Featurize mols in SDF, SMI.
    SMILES are regenerated with RDKit to get canonical SMILES without chirality information.

    Parameters
    ----------
    in_file : str
        Input SDF, SMI, ISM (or GZ)
    out_file : str
        Output csv
    model_path : str
        File path to pre-trained Gensim word2vec model
    r : int
        Radius of morgan fingerprint
    uncommon : str
        String to used to replace uncommon words/identifiers while training. Vector obtained for 'uncommon' will be used
        to encode new (unseen) identifiers

    Returns
    -------
    """
    # Load the model
    word2vec_model = word2vec.Word2Vec.load(model_path)
    if uncommon:
        try:
            word2vec_model[uncommon]
        except KeyError:
            raise KeyError('Selected word for uncommon: %s not in vocabulary' % uncommon)

    # File type detection
    in_split = in_file.split('.')
    f_type = in_split[-1].lower()
    if f_type not in ['sdf', 'smi', 'ism', 'gz']:
        raise ValueError('File extension not supported (sdf, smi, ism, sdf.gz, smi.gz)')
    if f_type == 'gz':
        if in_split[-2].lower() not in ['sdf', 'smi', 'ism']:
            raise ValueError('File extension not supported (sdf, smi, ism, sdf.gz, smi.gz)')
        else:
            f_type = in_split[-2].lower()

    print('Loading molecules.')
    if f_type == 'sdf':
        df = PandasTools.LoadSDF(in_file)
        print("Keeping only molecules that can be processed by RDKit.")
        df = df[df['ROMol'].notnull()]
        df['Smiles'] = df['ROMol'].map(Chem.MolToSmiles)
    else:
        df = pd.read_csv(in_file, index_col=False, usecols=[2, 0])[['Smiles', 'ID']]
        print(df.head(10))
        PandasTools.AddMoleculeColumnToFrame(df, smilesCol='Smiles')
        print("Keeping only molecules that can be processed by RDKit.")
        df = df[df['ROMol'].notnull()]
        df['Smiles'] = df['ROMol'].map(Chem.MolToSmiles)  # Recreate SMILES

    print('Featurizing molecules.')
    df['mol-sentence'] = df.apply(lambda x: MolSentence(mol2alt_sentence(x['ROMol'], r)), axis=1)
    vectors = sentences2vec(df['mol-sentence'], word2vec_model, unseen=uncommon)
    df_vec = pd.DataFrame(vectors, columns=['%03i' % x for x in range(vectors.shape[1])])
    df_vec.index = df.index
    df = df.join(df_vec)

    df.drop(['ROMol', 'mol-sentence'], axis=1).to_csv(out_file)