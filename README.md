## Usage

```shell
# initialize iLCSoft/key4hep
$ source /cvmfs/ilc.desy.de/key4hep/setup.sh

# compile this processor
$ mkdir build && cd build
$ cmake ..
$ make install -j20

# add processor to MARLIN_DLL
$ export MARLIN_DLL=$MARLIN_DLL:/path/to/source/lib/libPeterProcessor.so

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

NOTE: To run code without event display change eventDisplay parameter
in the xml steering file to false. No need to recompile.
