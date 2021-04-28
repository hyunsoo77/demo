# Python_plot

This repository is only for reviewing the example of code encapsulation for reproducibility.

Python_plot reads a tab-delimited file to generate a heatmap. This simple script will be encapsulated into a docker image. The docker image contains other dependent software including python3. Thus, python_plot.py inside this docker image can be executed without any software version issue.


### Prerequisites
You must have [Docker](https://www.docker.com/) installed in order to run python_plot.py with docker. No other dependencies are required. Please see the Docker documentation for installation information. If you don't have Docker, simply execute the following line after loading python3 environment (e.g. module load python/3.8.3):

```
python_plot.py --xlabel xx --ylabel yy --title title --figsize_y 5 heatmap input/table_mtx.txt output/table_mtx_heatmap.pdf
```

If you see any error message here, please make sure that you are using python3 environemnt and the dependent packages (os, sys, psutil, re, ssl, math, time, datetime, glob, ujson, shelve, numpy, matplotlib, seaborn, pandas, natsort, optparse, collections, itertools) were installed. Using an docker image significantly reduces the efforts to install these dependent packages; see [below](#installing).

The output file is located at output/table_mtx_heatmap.pdf.

```
├── output
│   └── table_mtx_heatmap.pdf
```

If you have installed acroread in your linux box, you can check the output PDF file. Otherwise, you may use your own PDF reader.

```
acroread output/table_mtx_heatmap.pdf
```


### Installing

The software is meant to be run as a Docker image. The docker image can be built by the following commands:
```
git clone https://github.com/hyunsoo77/demo.git
cd demo
docker build . -t hyunsoo77/hkim:demo
```

If you don't want to build docker image, you can download my image from DockerHub:
```
docker pull hyunsoo77/hkim:demo
```




## Command-line call

For example, if your current working directory ($PWD) looks like this:

```
.
├── Dockerfile
├── input
│   ├── table_2col_colors.txt
│   ├── table_2row_colors.txt
│   ├── table_col_colors.txt
│   ├── table_mtx.txt
│   └── table_row_colors.txt
├── output
│   └── table_mtx_heatmap.pdf
├── python_plot.py
└── README.md
```

where `outputs` is the directory you wish python_plot.py to place its results files, the command you would use to run python_plot.py inside docker image hyunsoo77/hkim:demo could be:

```
docker run -v $PWD/input:/home/hkim/in -v $PWD/output:/home/hkim/out hyunsoo77/hkim:demo --xlabel xx --ylabel yy --title title --figsize_y 5 heatmap in/table_mtx.txt out/table_mtx_heatmap.pdf
```


### File generated

The output file is a PDF file of heatmap generated by python_plot.py with $PWD/input/table_mtx.txt.


```
├── output
│   └── table_mtx_heatmap.pdf
```


