# SACO (Sequence Alignment COmpressor) #
SACO: a lossy compression tool for the sequences alignments found in the MAF files.

This compression tool was designed to handle the DNA bases and gap symbols that can be found in MAF files. Our method is based on a mixture of finite-context models. Contrarily to a recent approach ([Hanus 2010](http://dx.doi.org/10.1109/TIT.2009.2037052)), it addresses both the DNA bases and gap symbols at once, better exploring the existing correlations. For comparison with previous methods, our algorithm was tested in the [multiz28way](http://hgdownload-test.cse.ucsc.edu/goldenPath/hg18/multiz28way) dataset. On average, it attained 0.94 bits per symbol, approximately 7% better than the previous best, for a similar computational complexity. We also tested the model in the most recent dataset, [multiz46way](http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/multiz46way). In this dataset, that contains alignments of 46 different species, our compression model achieved an average of 0.72 bits per MSA block symbol.


# INSTALLATION #
In order to compile the source code, you will need to install a GCC compiler on a Unix platform (Linux or OS X). If you are using Windows, it will be easy to use the pre-compiled binaries that are in folders "bin/win32" and "bin/win64".

## Linux ##
For Linux users, install the build-essentials package which contains GCC and other utilities in order to be able to compile the source code. To install the build-essentials package type:
<pre>sudo apt-get install build-essential</pre>
After that you only need to type:
<pre>make -f Makefile.linux</pre>
to create the binaries "bin/linux/SACOe" (encoder) and "bin/linux/SACOd" (decoder).

## OS X ##
For OS X users, it depends on which Xcode version is installed. For the most recent versions, you will need to install the "Command Line Tool" in order to have the "make" utility. It seems that the "Command Line Tools" are not installed by default anymore when you install Xcode. In order to install them, open Xcode, go to Preferences -> Downloads -> Components -> Command Line Tools. This also should install a GCC compiler as well. If you want a recent compiler you can install it using Homebrew by typing the following command in a Terminal:
<pre>brew install gcc48</pre>
After that, we need to make sure that the "CC" variable in the "Makefile.osx" file is linked to the GCC previously installed. The most recent versions of XCode come with a modified version of GCC known as LLVM. This tool was not tested using LLVM so it will probably not work if you try to compile the code using it. In order to generate the binaries just type:
<pre>make -f Makefile.osx</pre>
to create the binaries "bin/osx/SACOe" (encoder) and "bin/osx/SACOd" (decoder).

## Windows ##
The source code was NOT tested in a Windows enviroment. Nevertheless, you can compile the code using a cross-compiler in a Linux environment after installing the cross-compiler [MinGW-w64](http://mingw-w64.sourceforge.net). After installing MinGW-w64, just type:
<pre>make -f Makefile.win32</pre>
to get the "bin/win32/SACOe32.exe" (encoder) "bin/win32/SACOd32.exe" (decoder) executables (32-bits architecture) and for the 64-bits architecture just type:
<pre>make -f Makefile.win64</pre> 
to get the "bin/win32/SACOe64.exe" (encoder) and "bin/win32/SACOd64.exe" (decoder) executables.
The encoder seems to work just fine however it seems that there is a bug in the decoder that will be fixed in the next commit...

# USAGE #

## Encoding ##
The SACOe, SACOe32.exe, and SACOe64.exe programs have several parameters that can be defined by the user. In the following you can find a description with the most relevant parameters available.

<pre>Usage: SACOe [options] ... [MAF File]</pre>
The most relevant options are:
<table align="center">
        <tr> 
          <td width="35%">-v</td> 
          <td width="65%">Activates vervose mode.</td>
        </tr>
        <tr> 
          <td width="35%">-h</td> 
          <td width="65%">Prints some help information.</td>
        </tr>
        <tr> 
          <td width="35%">-o [encodedFile]</td> 
          <td width="65%">If present, it writes the encoded data into file "encodedFile".</td>
        </tr>
        <tr> 
          <td width="35%">-e</td> 
          <td width="65%">Estimation only. Does not create the binary compressed file.</td>
        </tr>
        <tr> 
          <td width="35%">-alm</td> 
          <td width="65%">Activate the acenstral line mode.</td>
        </tr>
        <tr> 
          <td width="35%">-scm</td> 
          <td width="65%">Activate the static column model.</td>
        </tr>
        <tr> 
          <td width="35%">-cm1 [n/d t=threshold]</td> 
          <td width="65%">Columnwise Model 1.</td>
        </tr>
         <tr> 
          <td width="35%">-cmn [n/d t=threshold]</td> 
          <td width="65%">Columnwise Model N.</td>
        </tr>
         <tr> 
          <td width="35%">-u 0 [leftSize-rightSize n/d t=threshold]</td> 
          <td width="65%">Ancestral context model with "leftSize" symbols on the left and "rightSize" symbols on the right.</td>
        </tr>
        <tr> 
          <td width="35%">-u template [n/d t=threshold]</td> 
          <td width="65%">2D image context template. Templates available 1-14 and 20-24.</td>
        </tr>
        <tr> 
          <td width="35%">-g [gamma]</td> 
          <td width="65%">Gamma value used in the model mixture.</td>
        </tr>
  </table>

## Decoding ##

## Examples ##

# CITE #
If you use this software, please cite the following publications: 
* [Luís M. O. Matos](http://sweet.ua.pt/luismatos), [Diogo Pratas](http://sweet.ua.pt/pratas), and [Armando J. Pinho](http://sweet.ua.pt/ap), ["A Compression Model for DNA Multiple Sequence Alignment Blocks"](http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6415270), in [IEEE Transactions on Information Theory](http://ieeexplore.ieee.org/xpl/RecentIssue.jsp?punumber=18), volume 59, number 5, pages 3189-3198, May 2013.
* [Luís M. O. Matos](http://sweet.ua.pt/luismatos), [Diogo Pratas](http://sweet.ua.pt/pratas), and [Armando J. Pinho](http://sweet.ua.pt/ap), ["Compression of whole genome alignments using a mixture of finite-context models"](https://dl.dropboxusercontent.com/u/1944285/publications/ConferencePapers/Matos-2012c.pdf), in [Proceedings of the International Conference on Image Analysis and Recognition, ICIAR 2012](http://www.aimiconf.org/iciar12), Editors A. Campilho and M. Kamel, volume 2324 of Lecture Notes in Computer Science (LNCS), pages 359-366, Springer Berlin Heidelberg, Aveiro, Portugal, June 2012.

# ISSUES #
The windows decoders (SACOd32.exe and SACOd64.exe) have a bug that will be fixed very soon...
For other issues please use the [issues link](https://github.com/lumiratos/saco/issues) at GitHub.
<!---
At the time, there are no relevant issues detected but if you find one please let me know using the [issues link](https://github.com/lumiratos/saco/issues) at GitHub.
-->

# COPYRIGHT #
Copyright (c) 2014 Luís M. O. Matos. See [LICENSE.txt](https://github.com/lumiratos/saco/blob/master/LICENSE.txt) for further details.
