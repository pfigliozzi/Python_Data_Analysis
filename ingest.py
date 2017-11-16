import pandas as pd
import string
from sklearn.preprocessing import LabelEncoder, OneHotEncoder


class IngestFrame(pd.DataFrame):
    def __init__(self, data=None, index=None, columns=None, dtype=None, copy=False):
        pd.DataFrame.__init__(self, data, index, columns, dtype, copy)

        self._cat_to_label_encoders = None

    
    def encode_cats_to_labels(self):
        
        if self._cat_to_label_encoders is not None:
            raise Exception('Category encoders for training set already defined!')
            return
        self._cat_to_label_encoders = {}
        for column in self:
            if self[column].dtype == 'object':
                self.loc[pd.isnull(self[column]), column] = 'NaN'
                encoder = LabelEncoder()
                encoder.fit(self[column])
                self[column] = encoder.transform(self[column])
                self._cat_to_label_encoders[column] = encoder
    
    def transform_cats_test_set(self, DataFrame):

        for column, encoder in self._cat_to_label_encoders.iteritems():
            try:
                transformed = encoder.transform(DataFrame[column])
            except ValueError as e:
                error_string = str(e)
                labels = string.split(error_string, ': ')[-1]
                raise Warning("Test data column '"+column+"' has new labels: "+labels)




        
