Total least squares method: Matlab toolbox (Functions and Demos)


Authors: Petras, I., Bednarova, D., Skovranek, T. and Podlubny, I.
Technical University of Kosice, Faculty of BERG, Slovak Republic,
{ivo.petras, dagmar.bednarova, tomas.skovranek, igor.podlubny}@tuke.sk


*List of new functions:
- function [Err, P] = fit_2D_data(XData, YData, vizualization)
- function [Err, P] = fit_2D_data_SVD(XData, YData, vizualization)
- function [Err, P] = orm(XData, YData, order)
- function [F, Srez, Scel] = statindexes(XData, YData, a, b)
- function [Err, min_param] = numerFminS(fun, p, LBa, UBa, xdata, ydata)
- function [Err, min_param] = numerFminD(fun, p, LBa, UBa, xdata, ydata, FTime, Tvz)
- function [I] = corrindex(XData, YData, YDataM, par_number)


*List of already published functions:
- function [Err,N,P] = fit_3D_data(XData, YData, ZData, geometry, visualization, sod)
  url: http://www.mathworks.com/matlabcentral/fileexchange/12395, 2006.
- function [x,fval,exitflag,output] = fminsearchbnd(fun, x0, LB, UB, options, varargin)
  url: http://www.mathworks.com/matlabcentral/fileexchange/8277, 2006.

Above functions 'fit_3D_data' and 'fminsearchbnd' will be installed to the
computer via install supporting package 'requireFEXpackage(FEXSubmissionID)'
created by Igor Podlubny, 2011 (Matlab Central File Exchange ID: #31069).


*List of supporting functions:
- function [min_dist,CP] = dist_dsearch(points, M, graph)
- function [d] = distance_nD(Point1, Point2)
- function [dy] = dmodel2D(t,y,a)
- function [yM] = model(xm, a)
- requireFEXpackage(FEXSubmissionID)
  url: http://www.mathworks.com/matlabcentral/fileexchange/31069, 2011.


*Source (data) files:
- dataLRM.txt
- dataNRM.txt
- dataERM.txt
- dataHeater.txt
- switzerland.giu


*List of demos:
- demo3Dplane
- demoIdent
- demoLRM
- demoNRM


*Requirements:
- Matlab with Statistics Toolbox

*Note: type 'help FunctionName' for more information about every function.