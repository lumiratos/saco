# SACO (Sequence Alignment COmpressor) #
SACO: a lossy compression tool for the sequences alignments found in the MAF files.

This compression tool was designed to handle the DNA bases and gap symbols that can be found in MAF files. Our method is based on a mixture of finite-context models. Contrarily to a recent approach ([Hanus 2010](http://dx.doi.org/10.1109/TIT.2009.2037052)), it addresses both the DNA bases and gap symbols at once, better exploring the existing correlations. For comparison with previous methods, our algorithm was tested in the [multiz28way](http://hgdownload-test.cse.ucsc.edu/goldenPath/hg18/multiz28way) dataset. On average, it attained 0.94 bits per symbol, approximately 7% better than the previous best, for a similar computational complexity. We also tested the model in the most recent dataset, [multiz46way](http://hgdownload-test.cse.ucsc.edu/goldenPath/hg19/multiz46way). In this dataset, that contains alignments of 46 different species, our compression model achieved an average of 0.72 bits per MSA block symbol.


# INSTALLATION #

## Linux ##

## OS X ##

## Windows ##

# USAGE #

## Encoding ##

## Decoding ##

## Examples ##

# CITE #
If you use this software, please cite the following publications: 
* [Luís M. O. Matos](http://sweet.ua.pt/luismatos), [Diogo Pratas](http://sweet.ua.pt/pratas), and [Armando J. Pinho](http://sweet.ua.pt/ap), ["A Compression Model for DNA Multiple Sequence Alignment Blocks"](http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6415270), in [IEEE Transactions on Information Theory](http://ieeexplore.ieee.org/xpl/RecentIssue.jsp?punumber=18), volume 59, number 5, pages 3189-3198, May 2013.
* [Luís M. O. Matos](http://sweet.ua.pt/luismatos), [Diogo Pratas](http://sweet.ua.pt/pratas), and [Armando J. Pinho](http://sweet.ua.pt/ap), ["Compression of whole genome alignments using a mixture of finite-context models"](https://dl.dropboxusercontent.com/u/1944285/publications/ConferencePapers/Matos-2012c.pdf), in [Proceedings of the International Conference on Image Analysis and Recognition, ICIAR 2012](http://www.aimiconf.org/iciar12), Editors A. Campilho and M. Kamel, volume 2324 of Lecture Notes in Computer Science (LNCS), pages 359-366, Springer Berlin Heidelberg, Aveiro, Portugal, June 2012.

# ISSUES #
At the time, there are no relevant issues detected but if you find one please let me know using the [issues link](https://github.com/lumiratos/saco/issues) at GitHub.

# COPYRIGHT #
Copyright (c) 2014 Luís M. O. Matos. See [LICENSE.txt](https://github.com/lumiratos/saco/blob/master/LICENSE.txt) for further details.
