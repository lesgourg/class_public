import numpy as np 
#uses one dimensional interpolation
def PropagateDatawithListOld(k,FValue,zredindex,transferFunctionlist):
	return (transferFunctionlist[zredindex](k.ravel()) * FValue.ravel()).reshape(FValue.shape)

def PropagateDatawithList(k, FValue, zredindex, transferFunctionlist):
    result = {}
    for field, transfer_function in transferFunctionlist.items():
        result[field] = (transfer_function[zredindex](k.ravel()) * FValue.ravel()).reshape(FValue.shape)
    return result

#module with uses two dimensional interpolation and propagates all data at once (fastest but high memory consumption)
def PropagateAllData(k,FValue,allzred,transferFunction):

	allFValue = np.ones((len(allzred),FValue.shape[0],FValue.shape[1]),dtype=complex)

	for kxindex in range(FValue.shape[0]):
		allFValue[:,kxindex,:] = transferFunction(allzred,k[kxindex])*FValue[kxindex]


	return allFValue


#module with uses 2 dimensional interpolation (slowest but can be useful if the set of redshift changes very often)
def PropagateData(k,FValue,zred,transferFunction):

	FValuenew = np.ones(FValue.shape,dtype=complex)

	for kxindex in range(FValue.shape[0]):
		allFValue[kxindex,:] = transferFunction(zred,k[kxindex])*FValue[kxindex]


	return allFValue




