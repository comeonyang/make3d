

INSTALLATION:

* Unzipping lmirank.zip creates a directory \lmirank
* Place this directory in your MATLAB path
* If you would like to use lmirank without specifying initial
  conditions, either ensure SeDuMi is installed or else
  call lmirank via YALMIP. In the latter case, you still need 
  to have a SDP solver installed that is recognized by YALMIP.


TESTING:

* To test that lmirank.m is working, run lmiranktest.m:

    [At,c,K,yfeas,y,info] = lmiranktest; 

  lmiranktest will generate a feasible problem and try to 
  solve it using lmirank.
* If you have both SeDuMi and YALMIP installed, a second test 
  can be run:
	
    [A,Bperp,Ctperp,y,info,X,Y] = lmiranktest2;

  Type

    help lmiranktest2

  for an explanation of this test.




Feedback is very welcome.


Robert Orsi, Australian National University.
robert.orsi@anu.edu.au