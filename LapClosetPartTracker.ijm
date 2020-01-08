// Author: Sébastien Tosi (IRB Barcelona)
// Contact: sebastien.tosi@irbbarcelona.org
// Version: 1.0
// Date: 27/12/2019

// Path to input image and results
inputDir = "C:\\Users\\Seb\\Desktop\\in";
outputDir = "C:\\Users\\Seb\\Desktop\\out";

// Parameters for Cytopacq time-lapse
laprad = 9;
thr = 0.05;
maxlnkdst = 15;

// Read arguments from command line
arg = getArgument();
parts = split(arg, ",");
for(i=0; i<parts.length; i++) 
{
	nameAndValue = split(parts[i], "=");
	if (indexOf(nameAndValue[0], "input")>-1) inputDir=nameAndValue[1];
	if (indexOf(nameAndValue[0], "output")>-1) outputDir=nameAndValue[1];
	if (indexOf(nameAndValue[0], "laprad")>-1) laprad=nameAndValue[1];
	if (indexOf(nameAndValue[0], "thr")>-1) thr=nameAndValue[1];
	if (indexOf(nameAndValue[0], "maxlnkdst")>-1) maxlnkdst=nameAndValue[1];
}

images = getFileList(inputDir);
for(img=0; img<images.length; img++) 
{

setBatchMode(true);

// Open image
FileName = images[img];
open(inputDir+File.separator+FileName);

// Workflow
run("Properties...", "channels=1 slices=1 frames="+nSlices);
run("FeatureJ Laplacian", "compute smoothing="+d2s(laprad,0));
FilterID = getImageID();
for(i=0;i<nSlices;i++)
{
	selectImage(FilterID);
	setSlice(i+1);
	run("Find Maxima...", "noise="+d2s(thr,3)+" output=List exclude light");
	NewX = newArray(nResults);
	NewY = newArray(nResults);
	for(j=0;j<nResults;j++)
	{
		NewX[j] = getResult("X",j);
		NewY[j] = getResult("Y",j);
	}
	if(i==0)
	{
		newImage("Mask", "16-bit black", getWidth(), getHeight(), nSlices);	
		MaskID = getImageID();
		NObjs = lengthOf(NewX);
		for(j=0;j<lengthOf(NewX);j++)setPixel(NewX[j],NewY[j],j+1);
		LastX = NewX;
		LastY = NewY;
	}
	else
	{
		selectImage(MaskID);
		setSlice(i+1);
		LastX_buf = newArray(NObjs);
		LastY_buf = newArray(NObjs);
		Buf = newArray(lengthOf(NewX));
		for(k=0;k<lengthOf(NewX);k++)Buf[k] = 1/0;
		for(j=0;j<NObjs;j++)
		{
			MinDst2 = 1/0;
			Mink = -1;
			for(k=0;k<lengthOf(NewX);k++)
			{
				Dst2 = (pow(NewX[k]-LastX[j],2)+pow(NewY[k]-LastY[j],2));
				if((Dst2<MinDst2)&&(Dst2<Buf[k])&&(Dst2<maxlnkdst*maxlnkdst))
				{
					MinDst2 = Dst2;
					Buf[k] = Dst2;
					Mink = k; 
				}
			}
			if(Mink>-1)
			{
				setPixel(NewX[Mink],NewY[Mink],j+1);
				LastX_buf[j] = NewX[Mink];
				LastY_buf[j] = NewY[Mink];
			}
		}
		LastX = LastX_buf;
		LastY = LastY_buf;
		// Check which particles were lost and ensure they will not be further linked
		getHistogram(values, counts, 65536);
		for(k=0;k<NObjs;k++)if(counts[k+1]==0)LastX[k] = 1/0;
	}	
}

// Save label mask
selectImage(MaskID);
run("Properties...", "channels=1 slices=1 frames="+nSlices);
run("Bio-Formats Exporter", "save="+outputDir+File.separator+FileName+" compression=Uncompressed");
//save(outputDir+File.separator+FileName);

run("Close All");
setBatchMode("exit & display");

}
