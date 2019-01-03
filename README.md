# Setting Up

For the current version, it works in CMSSW_9_1_0_pre1. It probably works in other version as well, but this is what I've gotten it working with. Since this is made to work with nanoAOD files, this probably doesn't have to change for each CMSSW version since constant n-tuple type.

```sh
cmsrel CMSSW_9_1_0_pre3
cd CMSSW_9_1_0_pre1/src
cmsenv
git clone https://github.com/dteague/Analyzer
cd Analyzer
make -j 10
```

# Changes

+ First version, updates will be written here

# How to Use Code

## Particle Definitions:

### PartDet

All information about the particles is derived from the PartDet (Particle Details) files. An example of these files is found in src_sh{Analyses/example}, and the Analyses folder is made to put future analyses so it is easy to run over them. 

The code natively will run the analysis in the PartDet folder. To change this, run with the -C option and the location of these files.

### Info files

All of the files are in json format. Json has very particular formatting restrictions, so errors may come up because of this. It is suggested to put each json file through a json file checker if a json error crops up (may be added in future editions).

Types are infered by the contents, so it should be smart enough to figure out what you need. The basic particles are read in by the Particle class while the other files must be read in by the analyzer code.

| Filename              | file contents                           |
|-----------------------|:---------------------------------------:|
| Muon_info.json        | Muon cuts                               |
| Electron_info.json    | Electron cuts                           |
| FatJet_info.json      | Fat Jet (W jet) cuts                    |
| Gen_info.json         | Gen level info (don't change)           |
| Jet_info.json         | Jet Cuts                                |
| Tau_info.json         | Tau Cuts (no implimented)               |
|-----------------------|-----------------------------------------|
| Hist_entries.in       | Histograms in each folder               |
| Hist_syst_entries.in  | Histograms in each systematic folder    |
| Cuts.in               | Multiplicity cuts, also defines folders |
| Run_info.json         | Generic Run cuts (e.g. MET cuts)        |
| Systematics_info.json | systematics set up here                 |

Some files haven't been changed over to the json format and require a specific configuration.

The basic idea is particles read in their json files which have a "subparticle" name, eg, Muon1, Muon2 in the Muon_info.json file. In each subparticle set are the different cuts which can be defined. The cut structure first asks if a cut is going to be used, then the cut values are used. In the function src_C++{bset}, each cut that is used is put in a list; this is done to speed up implimentation/avoid redundances.

## Histograms Implimentations

Histgrams are defined in the Hist_entries.in file (as well as the syst version of this). The way histograms are read in is similar to how the actual ROOT objects take arguments, ie:

```
<NAME>  <BINS>  <MIN>  <MAX>   // OR
<NAME 2D>  <XBINS> <XMIN>  <XMAX>  <YBINS>  <YMIN>  <YMAX>
```

Each histgram made is put into each folder created. This may lead to redunances, but that's all we've got right now...

Many of the histograms can be grouped together, so to facilitate the removal process, blocks of similar histograms are grouped under a heading that starts with the keywork "Fill."  To remove the block, set the Fill heading to `false`.  Since the heading won't be seen by the program, the calculates done by the block won't be done either, so marginal speed gains will be made by program (less 100th of total time, so not signicant)

### Making a new Histogram

Adding a new histogram is fairly easy since the program is dynamic enough to hold most changes.  Two main things need to be done.

1. The histogram and information needs to be put into PartDet/Hist_info.in.  This includes name, bins, min, max as well as which heading the histogram will be stored under.  This can follow the template of the other histograms, so this is relatively easy
2. The histogram needs to be filled with the right values.  The filling of the histograms is done in the method `fillFolder`.  In this method, there are several if blocks for the different headings.  Go to the appropriate heading (or make one if a new one was made in the Hist_info.in file), and use the fill command, ie

```C++
histAddVal(<Value>, <Short name>)
```

The "Short name" is a genericized name to facilitate filling similar histograms. Its made by removing the Particle name, eg Muon1_Pt -> Pt. For multiparticle where you want a specific particle info, the generized name is turned into Part1 or Part2, eg for DiMuon pair, you want the leading muon's Pt, the generized name is Part1_Pt (note: the first particle is alwasy the leading particle).

## Folder Structure

Folders in the program are made when reading PartDet/Cuts.in.  By default, the program will always make the last significant cut (range is not [0,-1]) into a folder.  To add folders, simply put `***` before the cut without any space.

e.g.
```
NRecoMuon1               0  -1
NRecoTau1                2   2
***NDiTauCombinations    1   0
NSusyCombinations        1  -1
NDiJetCombinations       0  -1
```
In this example, there is a cut on Tau1, DiTaus, and a VBF cut.  The folders created are NDiTauCominations and NSusyCombinations (last significant cut).

The order of the cuts can also be rearranged if one wants to see cut flow in a different way.

# FAQ

+ Q: [What happens when program crashes?]
+ Q: [How Do I Add new "Selections?"]
+ Q: [How Do I Run the dang code?]
## What happens when program crashes?

Try to get more info:

```sh
make clean; DEBUG=1 make -j8
# Now run in gdb
gdb --args ./Analyzer -in <inputfile> -out output.root
run
# when it crashes
where
```

The where command should give you the function and inputs to the function where something went wrong. This may take more investigation, but this can save a lot of time identifying where the issue is taking place, and thus how to diagnose it

Here is an example of the gdb output:

```sh
The lines below might hint at the cause of the crash.
If they do not help you then please submit a bug report at
http://root.cern.ch/bugs. Please post the ENTIRE stack trace
from above as an attachment in addition to anything else
that might help us fixing this issue.
===========================================================
...
...
#12 Analyzer::getGoodRecoLeptons (this=this
entry=0x7fff69207790, lep=..., ePos=ePos
entry=CUTS::eRTau1, eGenPos=eGenPos
entry=CUTS::eGTau, stats=...) at src/Analyzer.cc:561
#13 0x0000000000463110 in Analyzer::preprocess (this=0x7fff69207790, event=0) at src/Analyzer.cc:130
#14 0x000000000041d4ac in main ()
===========================================================
```

As we can see, this error happened in src_sh{getGoodRecoLeptons}, specifically for the Tau1. While this isn't telling what the error is, it certainly can narrow down where the problem is happening (e.g. a tau related cut is set wrong)

## How Do I Add new "Selections?"

The code is set up in the following way:

+ FillCuts goes through the multiplicity cuts defined in the Cuts.in file. It checks the numbers from the file and if this event has particles within the range specified
+ getGoodParticles is run with run each respective getGood function and has an associated CUTS object. The getGood functions find the good Particles and put them into an array of the good particles.
+ fillFolder runs through all of the Fill groups to put the data into the histograms.

There are a few distinct parts that need to be addressed to get a new selection created.

### CUTS must be created

The actual cut is done on a value specified in the `Cuts.in` file. The naming convention is to call it "N<Cut Name>", but this doesn't have to be followed. Once a variable is created in the `Cuts.in` file, it needs to be linked to the code. This is done is the CUTS enum objects.

In the Cut_enum.h file, one just needs to put their enum into the list, before the value labelled Last. To link the CUTS variable with the `Cuts.in` variable, in `Analyzer.cc`, in the variable `cut_num`, you must create a map value in the vane of the others.

### getGood function must be created

Most likely, the new selection will require a unique way of applying cuts, so a new function must be created to accommodate this. One can use one the existing getGood functions as a template. The necessary part is simply that the function includes
```c
if(passCuts) active_part->at(ePos)->push_back(i);
```
or some variant. This line meaning the active_part (list of good parts being considered then), specifically the particles labelled by the CUTS tag of `ePos` is saying particle i is a good one, or it passed all of the cuts.

One will probably need to deal with the systematics which can be checked with the usual `needSyst(int)` function for Particle objects. If this doesn't make sense to you/you don't need systematics right now, just put at the top of the function:
```c
if(syst != 0) return;
```
So systematics aren't considered.

### A fillGroup must be created

After the above steps, the cut is made and applied. If one wants to plot specific variables connected to the cuts, one will have to create a fillGroup to handle this. 

First, histogram info is read from the PartDet file `Hist_entries.in`. It blocks related cuts into groups called filled groups, and to facilitate looking at these groups, groups of variables can be added or removed from the final root file by turning on or off the flag in front of the fillGroup name (e.g. "FillElectron1")

All of the histogram variables are put in the `fill_Folder` function in `Analyzer.cc`. There is a specific object for handling the different filling called FillVals which are structs in the `FillInfo.h` file (naming convention??). 

With this in mind, to get the histograms to work, one needs to create the cuts in the `Hist_entries.in` file. Then, a new mapped value needs to be created in `create_fillInfo()` in `Analyzer.cc`. The values put into the FillVals are used to make some things more generic (e.g. if you have `FILLER::Single`, it will all run through single particle filling). If more things needs to be added, change this `FillInfo.h`.

Last, in `fill_Folder`, there is a chain of if else blocks based on group name. Now, put filling into this if else and it should fill up the root file with the usual `histAddVal` function.

### In Summary

+ Put new cut name in `Cuts.in`
+ Create CUTS variable in `Cut_enum.h`
+ Link cut name and CUTS variable in `Analyzer.cc` in variable `cut_num`
+ create getGood function to fill up particle array at the right CUTS point
+ FillGroup created in `Hist_entries.in`
+ FillGroup linked to FillVals variable in `Analyzer.cc` in function `create_fillInfo()`
+ filling code put into if/else block in `Analyzer.cc` in function `fill_Folder`

## Q: How Do I Run the dang code?
The code works with options, so it should help you out if you need something done. Here is an example
```sh
./Analyzer -C Analyses/SSlep/ -in in_files/ -out test.root
```
+ The -C flag changes the directory where the particle information files are found. Normally, it defaults to ```PartDet```
+ The -in flag will read in a list of root files or read root files from a folder. In this example, a folder called `in_files/` was used
+ -out specifies the name of the outputted root file of the function.
