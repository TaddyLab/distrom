# distrom 

This is an R package for distributed multinomial regression.  The main function is `dmr', which fits independent Poisson log regressions for each category in your multinomial response.

To cite this package, use "Taddy (2013), Distributed Multinomial Regression".
The associated paper is available at <http://faculty.chicagobooth.edu/matt.taddy/research>.

The R package makes use of the 
    parallel
library.  This allows easy in-memory parallelization, and can also be used to facilitate distribution across multiple machines.

Any questions or comments can be emailed to <taddy@chicagobooth.edu>.

