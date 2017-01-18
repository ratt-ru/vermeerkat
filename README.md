# VerMeerKAT

RATT/RARG MeerKAT continuum self-calibration pipeline.

## Usage

### Download and image observations in the last 3 days

```bash
python -m vermeerat
```


### Specify a different configuration file

```bash
python -m vermeerat -c myconfig.conf
```

### Download and image a specific observation file, using a custom configuration file


```bash
python -m vermeerkat -f 123456789.h5 -c myconfig.conf
```

### The Astronomer, by Vermeer

<img src="https://upload.wikimedia.org/wikipedia/commons/0/0e/Johannes_Vermeer_-_The_Astronomer_-_WGA24685.jpg" alt="The Astronomer" align="middle" width="500" height="500"/>

### The latest version of the pipeline is depicted here. Unimplemented steps are shown in red:
![Pipeline](https://github.com/ska-sa/vermeerkat/blob/master/misc/Vermeerkat_flow.png)
