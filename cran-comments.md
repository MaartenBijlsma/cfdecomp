
# cfdecomp v0.4.0

This is the fourth release of the cfdecomp package. Compared to version 3, the cluster.resample function now makes use of rep() instead of a for() loop, speeding up the function considerably!

## Test environments (re-tested on August 11, 2021)
* local test on Windows 7 (64 bit), R version 4.1.0 (2021-05-18)


* Via check_win_release()
  * Ubuntu Linux 20.04.1 LTS, R-release, GCC
  * Fedora Linux, R-devel, clang, gfortran
  * Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* Via check_win_devel()
  * R Under development (unstable) (2021-08-09 r80724)
* via devtools::check_rhub()
  * Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  * Ubuntu Linux 20.04.1 LTS, R-release, GCC
  * Fedora Linux, R-devel, clang, gfortran



# cfdecomp v0.3.0

This is the third release of the cfdecomp package, which fixes a bug in the cluster.resample function.

## Test environments (re-tested on February 2, 2021)
* local test on Windows NT 6.1 (64-bit), R version 4.0.2 (2020-06-22)
* Via check_win_devel()
  * R Under development (unstable) (2021-01-31 r79912), 64 bit
* Via check_win_release()
  * R version 4.0.3 (2020-10-10), 64 bit
* via devtools::check_rhub()
  * Windows Server 2008 R2 SP1, R-devel, 32/64 bit
  * Ubuntu Linux 20.04.1 LTS, R-release, GCC
  * Fedora Linux, R-devel, clang, gfortran
* via R travis-ci Ubuntu 16.04.6 LTS
	* R version 4.0.2 (2020-06-22)

## R CMD check result
All of the above except one were OK and returned

0 errors | 0 warnings | 0 notes 

One (Ubuntu Linux 20.04.1 LTS, R-release, GCC) returned 0 errors, 0 warnings, and 1 note:

"Examples with CPU (user + system) or elapsed time > 5s"

The elapsed time was between 6 and 7 seconds.




# cfdecomp v0.2.0

This is the second release of the cfdecomp package, which adds two functions and some functionality to existing functions.

## Test environments (re-tested on April 30, 2020)

* local test on Windows NT 6.1 (64-bit), R version 3.6.2 (2019-12-12)
* Via check_win_devel()
  * Windows Server 2008 R2 SP1, R-devel (R Under development (unstable) (2020-04-30 (r78335)), 32/64 bit
* via devtools::check_rhub()
	* Ubuntu Linux 16.04 LTS, R-release, GCC
	* Fedora Linux, R-devel, clang, gfortran
* via R travis-ci Ubuntu 16.04.6 LTS
	* R version 4.0.0 (2020-04-24)

## R CMD check result
All of the above were OK and returned

0 errors | 0 warnings | 1 notes 

  New maintainer:
    Maarten Jacob Bijlsma <maarten.bijlsma@gmail.com>
  Old maintainer(s):
    Maarten Bijlsma <maarten.bijlsma@gmail.com>
    
This is in order: I normally include my middle name in publications. Not including it was an oversight in the first version. Since I have not yet advertised the package and it has not been referred to yet in publications, I hope this is OK.




# cfdecomp v0.1.0

This is a first release of the cfdecomp package, which has so far only been available on github under the name cfdecomp

## Comments received from CRAN team member at 1:57PM CET March 10, 2020 with response:

* DOI: specification is only supported in the DEWSCRIPTION file, otherwise you either have to use specific DOI markup such as \doi{} in Rd files or URL markup with fully specified URLs.
  * We now use the url in angle brackets in the README instead.


## Comments received from CRAN team member at 9:59AM CET on March 10, 2020 with response:

* Please shorten the title to a maximum of 65 characters. Acronyms/Abbreviations can be used on their own in the title as long as they are explained in the description field.
  * The title has now been changed to 61 characters (spaces included in the count).

* Please only capitalize names, sentence beginnings and
abbreviations/acronyms in the description text of your DESCRIPTION file.
F.i.: --> integration
  * 'Integration', 'Working' and 'Paper' are now no longer capitalized. We believe these were all words that were incorrectly capitalized.

* Please write references in the Description field of the DESCRIPTION file in the form authors (year) <doi:...>, authors (year) <arXiv:...>, authors (year, ISBN:...) or if those are not available: authors (year) <https:...>
with no space after 'doi:', 'arXiv:', 'https:' and angle brackets for auto-linking. (If you want to add the title as well, quote it. --> "title")
  * We now refer to Sudharsanan & Bijlsma (2019) with the digital object identifier added between angle brackets. Corresponding changes were made in README.md.

* Please always write TRUE and FALSE instead of T and F. (Also never name your variables T or F.)
  * All instances of =T were changed to =TRUE. Instances of =F were not found.

* \dontrun{} should be only used if the example really cannot be executed (e.g. because of missing additional software, missing API keys, ...) by the user. That's why wrapping examples in \dontrun{} adds the comment ("# Not run:") as a warning for the user. Please unwrap the examples if they are executable in < 5 sec, or create additionally small toy examples to allow automatic testing. You could also replace \dontrun{} with \donttest{}, but it would be preferable to have automatic checks for functions.
  * The examples have been unwrapped and altered to run on a subsample of the data and with fewer mc and bs iterations so that they are executable in < 5 seconds.

## Test environments (re-tested on March 10, 2020, after implementing changes described above)

* local test on Windows NT 6.1 (64-bit), R version 3.6.2 (2019-12-12)
* Via win-builder on x86_64-w64-mingw32 (64-bit)
  * R Under development (unstable) (2020-01-28 r77738)
* via devtools::rhub()
	* Windows Server 2008 R2 SP1, R-devel, 32/64 bit
	* Ubuntu Linux 16.04 LTS, R-release, GCC
	* Fedora Linux, R-devel, clang, gfortran
* via R travis-ci Ubuntu 16.04.6 LTS
	* R version 3.6.2 (2017-01-27)


## R CMD check result
All of the above were OK and returned

0 errors | 0 warnings | 1 notes 



NOTE: 
1) Maintainer: Maarten Bijlsma <maarten.bijlsma@gmail.com>
  
  New submission

This is in order.



2) Possibly mis-spelled words in DESCRIPTION:

    Bijlsma (11:671)
    Counterfactual (3:8)
    Sudharsanan (11:657)
    analytical (11:431)
    cfdecomp (11:14)
    counterfactual (11:231)
    multivariable (11:358) 

Bijlsma and Sudharsanan are the authors names
cfdecomp is the package name
Counterfactual and analytical are correctly written, as is multivariable. The latter is expressly NOT the same as 'multivariate'.


Many thanks

