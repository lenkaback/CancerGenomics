## Prerequisities using Python 3.8:
  - numpy 1.18.5
  - scikit-learn 0.23.2
  - pandas 1.1.2
  - matplotlib 3.3.2

These packages can be installed by `pip`:
```
pip3 install --user scikit-learn==0.23.2 pandas==1.1.2 numpy==1.18.5 matplotlib==3.3.2
```

## Results
The notebook `Deletions_ML_task.ipynb` can be run in this directory by using `jupyter notebook` or `jupyter lab`. The results are automatically written to this directory to `predictions_ada.tsv`. The model `ada_model.pickle` is saved here as well.

To load the model and predict on `data` (loaded as shown in the notebook), one can use:
```
import lzma
import pickle
with lzma.open("./ada_model.pickle", "rb") as model_file:
    ada = pickle.load(model_file)
predictions = ada.predict(data)
```
