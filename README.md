# Fluorescence Monte Carlo

This is a program designed for simulating fluorescence.

- In the GUI, you can set your parameters and click 'run', and wait for the result, the details will be printed in the command line.
- If you have more than one case to simulate, you can recored the parameters in a single csv file, then click 'File' > 'Auto' to load the file, and click 'Auto Run'.

## Dependencies
This program needs CUDA support, so make sure you have installed Nvidia Cuda environment.
If you want to open with GUI, you need to install python and run:

``` pip install -r requirements.txt ```

To build executable:

```make```


## Usage

To start GUI:

``` python FMC.py ```

![alt text](https://github.com/ToastCheng/temp/blob/master/fig1.png "Logo Title Text 1")

To use command line:
```./mcf [folder_name]```

note that you need to create a folder, say "example", and also there should be an txt file "example.txt" in the DRSOutput folder (the name should be the same). 

The format of example.txt should be:
wavelength    mua_top    mus_top    mua_bottom    mus_bottom

e.g.
```
400   2.90    150.01    1.34    130.23
410   3.56    140.33    1.56    124.23
420   1.21    120.22    1.00    101.46
...
```

If you use the GUI, then you don't have to worry about these things.

""
