function conversionmatrix = convertindexmatrix(formulation, indexmatrix)

rimodel = formulation.refractiveIndex.refractiveIndexModel;
rifun = rimodel.segments(2).fitting;

conversionRange = rimodel.segments(2).conversionRange;
conversion = (conversionRange(1):0.001:conversionRange(2))';
refractiveindex = feval(rifun,conversion);

conversionmatrix = interp1(refractiveindex,conversion,indexmatrix);