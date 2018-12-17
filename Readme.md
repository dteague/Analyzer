<div id="table-of-contents">
<h2>Table of Contents</h2>
<div id="text-table-of-contents">
<ul>
<li><a href="#sec-1">1. Setting Up</a></li>
<li><a href="#sec-2">2. Changes</a></li>
<li><a href="#sec-3">3. How to Use Code</a>
<ul>
<li><a href="#sec-3-1">3.1. Particle Definitions:</a>
<ul>
<li><a href="#sec-3-1-1">3.1.1. PartDet</a></li>
<li><a href="#sec-3-1-2">3.1.2. Info files</a></li>
</ul>
</li>
</ul>
</li>
<li><a href="#sec-4">4. FAQ</a></li>
</ul>
</div>
</div>

# Setting Up<a id="sec-1" name="sec-1"></a>

For the current version, it works in CMSSW<sub>9</sub><sub>1</sub><sub>0</sub><sub>pre1</sub>. It probably works in other version as well, but this is what I've gotten it working with. Since this is made to work with nanoAOD files, this probably doesn't have to change for each CMSSW version since constant n-tuple type.

    cmsrel CMSSW_9_1_0_pre1
    cd CMSSW_9_1_0_pre1/src
    cmsenv
    git clone https://github.com/dteague/Analyzer
    cd Analyzer
    make -j 10

# Changes<a id="sec-2" name="sec-2"></a>

-   First version, updates will be written here

# How to Use Code<a id="sec-3" name="sec-3"></a>

## Particle Definitions:<a id="sec-3-1" name="sec-3-1"></a>

### PartDet<a id="sec-3-1-1" name="sec-3-1-1"></a>

All information about the particles is derived from the PartDet (Particle Details) files. An example of these files is found in , and the Analyses folder is made to put future analyses so it is easy to run over them. 

The code natively will run the analysis in the PartDet folder. To change this, run with the -C option and the location of these files.

### Info files<a id="sec-3-1-2" name="sec-3-1-2"></a>

All of the files are in json format. Json has very particular formatting restrictions, so errors may come up because of this. It is suggested to put each json file through a json file checker if a json error crops up (may be added in future editions).

Types are infered by the contents, so it should be smart enough to figure out what you need. The basic particles are read in by the Particle class while the other files must be read in by the analyzer code.

<table border="2" cellspacing="0" cellpadding="6" rules="groups" frame="hsides">


<colgroup>
<col  class="left" />

<col  class="left" />
</colgroup>
<thead>
<tr>
<th scope="col" class="left">Filename</th>
<th scope="col" class="left">file contents</th>
</tr>
</thead>

<tbody>
<tr>
<td class="left">Muon<sub>info</sub>.json</td>
<td class="left">Muon cuts</td>
</tr>


<tr>
<td class="left">Electron<sub>info</sub>.json</td>
<td class="left">Electron cuts</td>
</tr>


<tr>
<td class="left">FatJet<sub>info</sub>.json</td>
<td class="left">Fat Jet (W jet) cuts</td>
</tr>


<tr>
<td class="left">Gen<sub>info</sub>.json</td>
<td class="left">Gen level info (don't change)</td>
</tr>


<tr>
<td class="left">Jet<sub>info</sub>.json</td>
<td class="left">Jet Cuts</td>
</tr>


<tr>
<td class="left">Tau<sub>info</sub>.json</td>
<td class="left">Tau Cuts (no implimented)</td>
</tr>
</tbody>

<tbody>
<tr>
<td class="left">Hist<sub>entries</sub>.in</td>
<td class="left">Histograms in each folder</td>
</tr>


<tr>
<td class="left">Hist<sub>syst</sub><sub>entries</sub>.in</td>
<td class="left">Histograms in each systematic folder</td>
</tr>


<tr>
<td class="left">Cuts.in</td>
<td class="left">Multiplicity cuts, also defines folders</td>
</tr>


<tr>
<td class="left">Run<sub>info</sub>.json</td>
<td class="left">Generic Run cuts (e.g. MET cuts)</td>
</tr>


<tr>
<td class="left">Systematics<sub>info</sub>.json</td>
<td class="left">systematics set up here</td>
</tr>
</tbody>
</table>

# FAQ<a id="sec-4" name="sec-4"></a>

-   Q: 3
-   Q:

-   A: <a id="What-happens-when-program-crashes" name="What-happens-when-program-crashes"></a>

Try to get more info:

    make clean; DEBUG=1 make -j8
    # Now run in gdb
    gdb --args ./Analyzer -in <inputfile> -out output.root
    run
    # when it crashes
    where

The where command should give you the function and inputs to the function where something went wrong. This may take more investigation, but this can save a lot of time identifying where the issue is taking place, and thus how to diagnose it

Here is an example of the gdb output:

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

As we can see, this error happened in , specifically for the Tau1. While this isn't telling what the error is, it certainly can narrow down where the problem is happening (e.g. a tau related cut is set wrong)

Folders

Folders in the program are made when reading PartDet/Cuts.in.  By default, the program will always make the last significant cut (range is not [0,-1]) into a folder.  To add folders, simply put \`\`\`\*\*\*\`\`\` before the cut without any space.

e.g.
\`\`\`
NRecoMuon1               0  -1
NRecoTau1                2   2
\*\*\*NDiTauCombinations    1   0
NSusyCombinations        1  -1
NDiJetCombinations       0  -1
\`\`\`
In this example, there is a cut on Tau1, DiTaus, and a VBF cut.  The folders created are NDiTauCominations and NSusyCombinations (last significant cut).

The order of the cuts can also be rearranged if one wants to see cut flow in a different way.

\### Histogram Management

All of the Histograms are stored in PartDet/Hist<sub>info</sub>.in.  On each line, the details of the histogram are:
\`\`\`
<NAME>  <BINS>  <MIN>  <MAX>   // OR
<NAME 2D>  <XBINS> <XMIN>  <XMAX>  <YBINS>  <YMIN>  <YMAX>
\`\`\`
Since the histogram information is read at the beginning of each run, the binning and domain of the histogram can be changed to fit the analysis.

As with all of the info files, the file supports C and python style line commenting (// and #).  This means, to remove a specific histogram, simply comment it out

Many of the histograms can be grouped together, so to facilitate the removal process, blocks of similar histograms are grouped under a heading that starts with the keywork "Fill."  To remove the block, set the Fill heading to 0 or false.  Since the heading won't be seen by the program, the calculates done by the block won't be done either, so marginal speed gains will be made by program (less 100th of total time, so not signicant)

\### New Histograms

Adding a new histogram is fairly easy since the program is dynamic enough to hold most changes.  Two main things need to be done.

1.  The histogram and information needs to be put into PartDet/Hist<sub>info</sub>.in.  This includes name, bins, min, max as well as which heading the histogram will be stored under.  This can follow the template of the other histograms, so this is relatively easy
2.  The histogram needs to be filled with the right values.  The filling of the histograms is done in the method \`\`\`fillFolder\`\`\`.  In this method, there are several if blocks for the different headings.  Go to the appropriate heading (or make one if a new one was made in the Hist<sub>info</sub>.in file), and write the following command to write to the histogram:
