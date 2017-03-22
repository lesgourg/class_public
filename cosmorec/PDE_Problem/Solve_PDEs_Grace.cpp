//==================================================================================================
// support for grace
//==================================================================================================
#include <grace_np.h>

//==================================================================================================
// grace output
//==================================================================================================
void prepare_Grace(double xmin, vector<double> &resonances, bool write_to_disk=0)
{
	if (GraceOpen(2048) == -1) 
	{
		cerr << " can't run grace " << endl; 
		exit(-1);
	}
	
    GracePrintf("PAGE SIZE 700, 500");
    
    //===============================================================
	// Set output device 
	//===============================================================
    if(write_to_disk)
    {
        //GracePrintf("PAGE SIZE 1500.0, 1250.0");
        GracePrintf("HARDCOPY DEVICE \"JPEG\"");
        //GracePrintf("DEVICE \"JPEG\" FONT on");
        GracePrintf("DEVICE \"JPEG\" FONT ANTIALIASING on");
        GracePrintf("DEVICE \"JPEG\" DPI 300");
        GracePrintf("DEVICE \"JPEG\" OP \"quality:100\"");
        //GracePrintf("DEVICE \"JPEG\" OP \"dct:float\""); 
        //GracePrintf("DEVICE \"JPEG\" OP \"optimize:on\""); 
        //GracePrintf("DEVICE \"JPEG\" OP \"smoothing:10\""); 
    }

	GracePrintf("XAXIS LABEL \"\\xn / n\\s\\021\"");
	GracePrintf("YAXIS LABEL \"\\1x\\S\\03\\N \\xD\\1n\\s\\1x\"");
	
	//===============================================================
	// Send some initialization commands to Grace 
	//===============================================================
	GracePrintf("world xmax 1.33");
	GracePrintf("world xmin 0.4");
	GracePrintf("world ymax 300");
	GracePrintf("world ymin 0.001");
	GracePrintf("XAXES SCALE NORMAL");
	GracePrintf("YAXES SCALE LOGARITHMIC");
	GracePrintf("TITLE SIZE 1");
	
	GracePrintf("xaxis tick major 0.1");
	GracePrintf("xaxis tick minor 0.05");
	GracePrintf("yaxis tick major 10");
	GracePrintf("yaxis tick minor 1");
	//
	GracePrintf("s0 LINEWIDTH 2");
	GracePrintf("s0 on");
	//
	GracePrintf("s2 COLOR 1");
	GracePrintf("s2 LINESTYLE 2");
	GracePrintf("s2 on");
	//
	
	GracePrintf("g0.s2 point %10.16f, %10.16f", -100.0, 1.0e-40);
	
	for(int m=0; m<(int)resonances.size(); m++)
	{
		if(m%2)
		{	
			GracePrintf("g0.s2 point %10.16f, %10.16f", resonances[m]*(1.0-1.0e-10), 1.0e-40);
			GracePrintf("g0.s2 point %10.16f, %10.16f", resonances[m]*(1.0+1.0e-10), 1.0e+40);
		}
		else
		{	
			GracePrintf("g0.s2 point %10.16f, %10.16f", resonances[m]*(1.0-1.0e-10), 1.0e+40);
			GracePrintf("g0.s2 point %10.16f, %10.16f", resonances[m]*(1.0+1.0e-10), 1.0e-40);
		}
	}
	
	return;
}


//==================================================================================================
// output figure
//==================================================================================================
void output_Grace(double zout, bool write_to_disk, vector<double> &xarr, 
                  vector<double> &yarr, string exten, int output_count)
{
	GracePrintf("KILL G0.S0");
	GracePrintf("KILL G0.S1");
	
	for (int i = 0; i < (int)xarr.size(); i++) 
	{
		double xxx=xarr[i];
		double yyy=yarr[i]*pow(xxx, 3);
		GracePrintf("g0.s0 point %10.16f, %10.16f", xxx, yyy);
	}
	
    string command;

    if(write_to_disk) 
    {
        string namef="./temp/JPEG/sol"+exten+"."+int_to_string(output_count, 5)+".jpg";
        command="PRINT TO \""+namef+"\"";
        GracePrintf(command.c_str());
    }
	
	std::ostringstream strs;
	strs.width(5);
	strs << zout;
	string dumstr= strs.str();
	
	command="TITLE \"HI Lyman-series spectral distortion at \\1z\\0 = "+dumstr+"\"";
	GracePrintf(command.c_str());

	if(write_to_disk) GracePrintf("PRINT");
	GracePrintf("redraw");
		
//	wait_f_r();

	return;
}			

//==================================================================================================
//==================================================================================================

