# 3D expanding pipe.

To build a 3D pipe with sudden expansion, do the following:

* Make certain you have built the nek tools via:

`maketools all`

* After making the tools (specifically, genmap, genbox, pretex), return to this working directory and type

`mkmesh`

* Then type

`makenek expansion`
`nek expansion`      # run on one core

or

`nekmpi expansion 4`      # run on 4 cores

or

`nekbmpi expansion 4`     # run on 4 cores in background mode

Of course, you will need nek, nekmpi, nekbmpi in your path.

You can change the SIZE file and recompile in the usual way.
If you increase lx1,etc., you probably should reduce lelt and
use more processors.

With increased resolution you probably also need to reduce dt,
which can be done by editing the base cyl2d.rea file prior to
mesh construction or by editing the resultant "expansion.rea" file.


The mesh construction procedure is described in the mkmesh script
in this directory.

By studying this script you should be able to see how to do basic
mesh refinement.   Many more options open up if you actually use the
graphical preprocessor, prenek, to edit the base file cyl2d.rea.
