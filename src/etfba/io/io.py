'''Difine IO functions.'''


from os.path import splitext
import re
from pickle import dump, load
from math import log, exp
import pandas as pd


def read_values(source, log_transform = False):
    '''
    Parameters
    ----------
    source: dict, Series or str
        If str, 'source' refers to the file name ending with .xlsx, .
        tsv, or .bin (pickled).
    log_transform: bool
        Indicates whether to perform natural logarithmic 
        transformation.
        
    Returns
    -------
    data: dict
    '''
    
    if isinstance(source, dict):
        data = source
    
    elif isinstance(source, pd.Series):
        data = source.to_dict()
    
    elif isinstance(source, str):
        ext = splitext(source)[1]
        if re.search(r'tsv', ext):
            data = pd.read_csv(
                source, 
                sep = '\t', 
                comment = '#', 
                header = None, 
                index_col = 0
            ).squeeze().to_dict()
        
        elif re.search(r'xls', ext):
            data = pd.read_excel(
                source, 
                comment = '#', 
                header = None, 
                index_col = 0,
            ).squeeze().to_dict()
        
        elif re.search(r'bin', ext):
            data = load(open(source, 'rb'))
    
    else:
        raise TypeError('can only read from dict, Series or file')
        
    if log_transform:
        data = {key: log(value) for key, value in data.items()}

    return data
    

def load_model(file):
    '''
    Parameters
    ----------
    file: str
        Filename ending with .bin.
    '''    
    
    return load(open(file, 'rb'))


def save_values(file, data, exp_transform = False):
    '''
    Parameters
    ----------
    file: str
        Filename ending with .xlsx, .tsv, or .bin to save the data.
    data: dict
        Data to be saved.
    log_transform: bool
        Indicates whether to perform natural exponential 
        transformation on the data.
    '''
    
    if exp_transform:
        data = {key: exp(value) for key, value in data.items()}
    
    ext = splitext(file)[1]
    if re.search(r'tsv', ext):
        if isinstance(list(data.values())[0], list):
            data = pd.DataFrame(data).T
        else:
            data = pd.Series(data) 
        data.to_csv(file, sep = '\t', header = False, index = True)
    
    elif re.search(r'xls', ext):
        if isinstance(list(data.values())[0], list):
            data = pd.DataFrame(data).T
        else:
            data = pd.Series(data)
        data.to_excel(file, header = False, index = True)
    
    elif re.search(r'bin', ext):
        dump(data, open(file, 'wb'))
    
    else:
        raise TypeError('can only save to .tsv, .bin or excel file')
    

def save_model(file, model):
    '''
    Parameters
    ----------
    file: str
        Filename ending with .bin to save the model.
    model: Model
        Model to be saved.
    '''

    dump(model, open(file, 'wb'))
