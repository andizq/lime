# Dust continuum and line emission for an envelope+disc model

Type the following commands:

```bash
$ python make_model.py
$ curl https://home.strw.leidenuniv.nl/~moldata/datafiles/ch3cn.dat -o ch3cn.dat
$ lime -nS -p 4 model.c   #Normal output (-n); sf3dmodels (-S) and 4 threads (-p 4)
```