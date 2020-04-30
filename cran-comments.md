
# cfdecomp v0.2.0

## Test environments (re-tested on April 30, 2020)

* local test on Windows NT 6.1 (64-bit), R version 3.6.2 (2019-12-12)
* Via win-builder on x86_64-w64-mingw32 (64-bit)
  * R Under development (unstable) (2020-01-28 r77738)
* via devtools::rhub()
	* Windows Server 2008 R2 SP1, R-devel, 32/64 bit
	* Ubuntu Linux 16.04 LTS, R-release, GCC
	* Fedora Linux, R-devel, clang, gfortran
* via R travis-ci Ubuntu 16.04.6 LTS
	* R version 3.6.2 (2017-01-27)




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

