## Usage

```shell
# initialize iLCSoft/key4hep
$ source /cvmfs/ilc.desy.de/key4hep/setup.sh

# compile this processor
$ mkdir build && cd build
$ cmake ..
$ make install -j20

# add processor to MARLIN_DLL
$ export MARLIN_DLL=$MARLIN_DLL:/path/to/source/lib/libTOFInfo.so

# run
$ Marlin ../xml/steer.xml

```

## Event display

Open a terminal on the local PC and do:  
```shell
$ source /afs/desy.de/project/ilcsoft/sw/x86_64_gcc75_ub1804/v02-02-02/init_ilcsoft.sh
$ glced -trust naf-ilc
```

Before running `Marlin ../xml/steer.xml` on the naf machine do:  
```shell
$ export CED_HOST=your_local_machine_name
# e.g.
# export CED_HOST=flc19
```

I have added this for debugging purposes, to clearly see a track/hits.
But the code will work without these steps just fine.


## What is inside

It shows how to extract ECAL hits and track information you will need for training.

In this example I just print it in terminal, but you can dump it into any other format you like.

The most important code is inside:

`printOutputFCN()`  
and  
`printOutputGNN()`

they are executed inside `processEvent()` which is the main function which executes every event.

## Output

After executing, this is the output example for `printOutputFCN()` with "Frank hits".

Note: you will need to scroll up a bit through the event display output.

```
[ VERBOSE "peterExample"] ***********************************************
[ VERBOSE "peterExample"] ****************Event 1**************************
[ VERBOSE "peterExample"] ***********************************************
[ VERBOSE "peterExample"] 
[ VERBOSE "peterExample"]         ======PFO 2===== PDG: 321 ======
[ VERBOSE "peterExample"] Track p: 3.84    pT: 2.04    pz: 3.253 (GeV)    theta: 32.1 (deg)    "TRUE" TOF: 9.579 ns
[ VERBOSE "peterExample"] Hit #    time true (ns)    time 50ps (ns)    energy (GeV)    d (mm)    d|| (mm)    dT (mm)
[ VERBOSE "peterExample"]     1         9.605          9.536             0.008542       7.733       7.599       1.434
[ VERBOSE "peterExample"]     2         9.617          9.516             0.009633       13.1       12.94       2.08
[ VERBOSE "peterExample"]     3         9.654          9.616             0.006505       23.58       23.54       1.286
[ VERBOSE "peterExample"]     4         9.666          9.667             0.006727       26.19       26.19       0.4609
[ VERBOSE "peterExample"]     5         9.704          9.701             0.00821       36.81       36.8       1.107
[ VERBOSE "peterExample"]     6         9.716          9.67             0.009072       42.17       42.13       1.743
[ VERBOSE "peterExample"]     7         9.753          9.784             0.008343       52.76       52.74       1.568
[ VERBOSE "peterExample"]     8         9.766          9.805             0.007452       55.41       55.39       1.516
[ VERBOSE "peterExample"]     9         9.803          9.753             0.00968       66.04       66       2.325
[ VERBOSE "peterExample"]     10         9.816          9.852             0.009868       71.38       71.33       2.6

```
