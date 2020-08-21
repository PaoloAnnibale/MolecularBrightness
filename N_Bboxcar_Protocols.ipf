//////////////////////////////////////////
/////////////////Readme///////////////////
//////////////////////////////////////////

//This procedure file contains the multiple functions necessary to analyze 
//fluorescence microscopy images and extract molecular number and brightness information from
//time series (temporal brightness) or still snapshots (spatial brightness)
//
//The functions actually dealing with brightness calculation are
//1. temporal_brighntess() which calulates brightness values from each pixel of a time series. The user has the option to select
//a polygonal ROI.This function can accomodate series of images acquired both with photon counting detectors (recommended) or
//analog detectors. All data treatment operations consist in the calculation of mean and variances.
//Upon user selection (doboxcar=1) the function can call
//2. boxcar2(), a function implementing a boxcar filter (as described in Trullo et al. 2013, 
//DOI 10.1002/jemt.22277) to remove the contribution of slow (relative to the boxcar window) 
//fluctiations to the variance of the time series, and hence the brightness. 
//The function contains also the code to provide a Gaussianity score for the selected ROI. This score is determined as the ratio
//of the actual score over the threshold (As defined upon performing the Kolmogorov-Smirnov (KS) goodness-of-fit test for two continuous distributions
//(that of actual pixels in the ROI and that of the pixel of an 'ideal' ROI, Gauss distributed with the same mean and variance as the first.)
//3. SpIDA_photoncounting() Calculates one brightness value from a single image (or ROI thereof) based on the relationship between the variance
//and the mean of the pixel intensities in the region. The code works for images acquired using photon counting detectors (recommended). For the analog
//case we refer here to Godin et al. PNAS 2011, doi/10.1073/pnas.1018658108 and to the related code available under https://neurophotonics.ca/software
//
//THe concepts contained in the above three functions can be exported to any other language/code, and provide the simple foundation of brightness calculation.
//
//The other functions, under Auxiliary Functions simply deal with the polygonal selection of the ROIs and with the graphical display of the results of the analysis, and are specific to IgorPro programming language. 
//
//Paolo Annibale, 2020


#pragma rtGlobals=1		// Use modern global access method.
#pragma rtGlobals=1		// Use modern global access method.
#include  <TransformAxis1.2>
#pragma rtGlobals=1		// Use modern global access method.
#include <Math Utility Functions>
//////////////////////////////////
///Creates acquisition Menu///////
//////////////////////////////////

Menu "Brightness menu"
	"temporal brightness", temporal_brighntess()
	"spatial brightness", SpIDA_photoncounting()
	"highligh on brightness plot", highlightBch1_poly()
	"highligh on intensity plot", highlightIch1_poly()
End

///////////////////////////////


///////////////////////////////////////////////////////
////////////calculates Temporal Brightness/////////////
///////////////////////////////////////////////////////
function temporal_brighntess()

variable startframe=0 //determine initial and final frames for the acquisition
variable endframe=0 //if last frame is put to zero, all frames are analyzed

variable doboxcar=1 //if 1, performs boxcar detrend according to Trullo et al. 2013

//these values characterize the detector performance.
variable offset=1E-7  //detector koff, i.e. dark counts baseline. Effectively 0 for photon counting
variable S=1 // Photon conversion factor for analog detectors. For Photon counting set to 1
variable sigma0=0 //width of the dark counts distribution. For photon counting set to 1
variable channel=1 //1 by default


		////////////////////////////
		////////GUI function////////
		////////////////////////////
		//prompts the user in order to read the string name
		string waveNM 
		
			
			Prompt waveNM,"current cell",popup,WaveList("*",";","DIMS:3")	// Set prompt for y param
			DoPrompt "Select current  cell name", waveNM
			if (V_Flag)
				return -1					// User canceled
			endif
			///////
		string/G wn=waveNM
		//////////

variable/G sizex, sizey, length
wave channel1=$(wn)
sizex=DimSize(channel1, 0)
sizey=DimSize(channel1, 1)
length=DimSize(channel1, 2)

	if(startframe!=0 && startframe>length)
		startframe=0
	
	endif

	if(endframe!=0)
		length=endframe
	
	endif

variable f, i, j, boxsize

boxsize=21 // length of boxcar filter window. Can be adjusted depending on experimentl settings (e.g. bleaching rate, motion of the cell...). Too short a window will suppress bona-fide molecular fluctuations


//initializes auxiliary vectors used throughout the calculations
make/O/N=(sizex, sizey) intensity, M_ROIMask
make/O/N=(sizex, sizey) sigma
make/O/N=(sizex, sizey) brightness //B dei papers
//make/O/N=(length-startframe) slice
make/O/N=(length-startframe) av_intensity


/////////////////////////////////////////////////////////////
////Creates and displays the average image from the stack////
/////////////////////////////////////////////////////////////
		ImageTransform/O/Meth=1 zProjection channel1
		dowindow/K average_intensity_preview
		display/N=average_intensity_preview /M/W=(30,25,45,35 ); appendimage/W=average_intensity_preview M_zProjection
		ModifyImage M_zProjection ctab= {*,*,Grays,0}
		ModifyGraph fSize=16,axThick=2

		////////////////////////////
		////////GUI function////////
		////////////////////////////
		
		////Prompts the user to select a region of interest where to perform the TB analysis
		
					DoWindow/F average_intensity_preview			// Bring graph to front
					
					if (V_Flag == 0)									// Verify that graph exists
						Abort "UserCursorAdjust: No such graph."
						return -1
					endif
				
					NewPanel/K=2 /W=(139,341,382,432) as "Pause for Cursor"
					DoWindow/C tmp_PauseforCursor_preview								// Set to an unlikely name
					AutoPositionWindow/E/M=1/R=average_intensity_preview		// Put panel near the graph
				
					DrawText 21,20,"Adjust the polygon cursors and then. Doubleclick when finished"
					DrawText 21,40,"Click Continue."
					
					Button button1,pos={20,58},size={92,20}, title="Draw Polygon "
					Button button1, proc=puntinelpoligono_preview
				
					Button button0,pos={120,58},size={92,20}, title="Continue"
					Button button0, proc=autochiuditi_preview
					
					PauseForUser tmp_PauseforCursor_preview, average_intensity_preview
		
					ImageBoundaryToMask width=sizex, height=sizey, xwave=x_poly, ywave=y_poly, seedX=mean(x_poly), seedY=mean(y_poly) 


/////////////////////////////////////////////////////////////
////Calculates the brightness image from the points//////////
////selected in the first ROI definition/////////////////////
/////////////////////////////////////////////////////////////

	for(i=0; i<sizex; i+=1)
	
		for(j=0; j<sizey; j+=1)
		
			if(M_ROIMask[i][j]==1)		

			ImageTransform/Beam={(i),(j)} getbeam channel1
			deletepoints 0, startframe, W_beam


				if(doboxcar==1)//calculates brightness using Boxcar function
					wave output=boxcar2(W_beam, boxsize, offset, S, sigma0) //this function (see below) does the actual brightness calculation from the boxcar windows
					intensity[i][j]=output[0]
					sigma[i][j]=output[1]
					brightness[i][j]=output[2]
				else //calculates brightness directly
					wavestats/Q W_beam
					intensity[i][j]=V_avg
					sigma[i][j]=V_sdev
					brightness[i][j]=(V_sdev^2-sigma0^2)/(V_avg-offset)
				endif
			else
					intensity[i][j]=0
					sigma[i][j]=1
					brightness[i][j]=1
			endif	
				
		endfor
	
	endfor

//plots the graphs and updates upon recall of the function and changes of S, sigma0, and offset

		make/O/N=(sizex, sizey) intensity_cal1, brightness_cal1, sigma_cal1
		intensity_cal1=(intensity-offset)/S
		brightness_cal1=brightness/S
		
		dowindow/K intensity_ch1
		display/N=intensity_ch1 /M/W=(30,25,45,35 ); appendimage/W=intensity_ch1 intensity_cal1
		ModifyImage intensity_cal1 ctab= {*,*,Grays,0}
		ModifyGraph fSize=16,axThick=2
		
		dowindow/K brightness_ch1 
		display/N=brightness_ch1/M/W=(0,0,15,10 ); appendimage/W=brightness_ch1 brightness_cal1
		ModifyImage brightness_cal1 ctab= {*,*,Spectrum,1}
		ModifyGraph fSize=16,axThick=2
		
		 
		dowindow/K BfracS_vs_I_ch1
		display/N=BfracS_vs_I_ch1/M/W=(0,15,15,25 ); appendtograph/W=BfracS_vs_I_ch1 brightness_cal1 vs intensity_cal1
		ModifyGraph mode=2
		ModifyGraph rgb=(65535,32768,32768)
		ModifyGraph fSize=16,axThick=2
		Label left "brightness"
		Label bottom "intensity"



end

///////////////////////////////////////////////////////////
/////////////////////SpIDA photon counting ////////////////
///////////////////////////////////////////////////////////

function SpIDA_photoncounting()


//detector calibration parameters
variable offset=0 //detector koff, i.e. dark counts baseline. Effectively 0 for photon counting
variable S=1 //Photon conversion factor for analog detectors. For Photon counting set to 1
variable sigma0=1 //width of the dark counts distribution. For photon counting set to 1
variable doscoring=1 //determines whether scoring for ROI homogeneity shall be performed.

//prompts the user in order to read the string name of SpIDA file
string waveNM 
	
	Prompt waveNM,"current cell",popup,WaveList("*",";","DIMS:2")	// Set prompt for y param
	DoPrompt "Select current  cell name", waveNM
	if (V_Flag)
		return -1					// User canceled
	endif
	///////
string/G wn_SpIDA=waveNM
//////////

variable i,j,k
variable/G sizex, sizey
wave SpIDA_wave=$(wn_SpIDA)
sizex=DimSize(SpIDA_wave, 0)
sizey=DimSize(SpIDA_wave, 1)

make/O/N=(sizex, sizey) M_ROIMask
make/o/N=(sizex*sizey) pixelsinROI

print WaveDims(Spida_wave)
	////Creates and display the SpIDA picture
		dowindow/K SpIDA_preview
		display/N=SpIDA_preview /M/W=(30,25,45,35 ); appendimage/W=SpIDA_preview SpIDA_wave
		//ModifyImage SpIDA_wave ctab= {*,*,Grays,0}

///////////////
////Prompts the user to select a region of interest where to perform the N&B analysis

			DoWindow/F SpIDA_preview		// Bring graph to front
			
			if (V_Flag == 0)									// Verify that graph exists
				Abort "UserCursorAdjust: No such graph."
				return -1
			endif
		
			NewPanel/K=2 /W=(139,341,382,432) as "Pause for Cursor"
			DoWindow/C tmp_PauseforCursor_preview			// Set to an unlikely name
			AutoPositionWindow/E/M=1/R=SpIDA_preview		// Put panel near the graph
		
			DrawText 21,20,"Adjust the polygon cursors and then. Doubleclick when finished"
			DrawText 21,40,"Click Continue."
			
			Button button1,pos={20,58},size={92,20}, title="Draw Polygon "
			Button button1, proc=puntinelpoligono_SpIDA_preview
		
			Button button0,pos={120,58},size={92,20}, title="Continue"
			Button button0, proc=autochiuditi_SpIDA_preview
			
			PauseForUser tmp_PauseforCursor_preview, SpIDA_preview

			ImageBoundaryToMask width=sizex, height=sizey, xwave=x_poly, ywave=y_poly, seedX=mean(x_poly), seedY=mean(y_poly) 
		
			k=0
			
			for(i=0; i<sizex; i+=1)
	
				for(j=0; j<sizey; j+=1)
		
					if(M_ROIMask[i][j]==1)	
					pixelsinROI[k]=SpIDA_wave[i][j]
					k+=1	
					endif				
				endfor
			
			endfor
			deletepoints (k), (sizex*sizey-k), pixelsinROI

			wavestats/Q pixelsinROI
			variable brightness=(V_sdev)^2/V_avg-1
			variable number=V_avg/brightness
			print "the brightness is: "+num2str(brightness)
			print "the number is: "+num2str(number)


			//this part scores each ROI for Gaussianity. Iterative ROI selection should maximize this value
			
		    	if(doscoring==1)
		    	make/O/N=(V_npnts) theory_data
		    	wave W_KSResults
		    	variable score
			 	theory_data=gnoise(V_sdev/sqrt(2)) + V_avg //generates Gaussian distributed intensity values for a ROI containing as many pixels as pixelsinROI and its mean and variance
		 		StatsKSTest/ALPH=0.05/Z/Q/T=1 pixelsinROI,theory_data //Calculates for Gaussianity using Komogorov-Smirnov criterion
	 	 		print "the gaussianity score is: "+num2str(W_KSResults[5]/W_KSResults[4]) //the bigger this value the better (it is the ratio between threshold and measured value, the measured value should be below threshold for gaussianity
			 	endif

end

///////////////////////////////////////////
/////////////Boxcar for images/////////////
///////////////////////////////////////////

function/wave boxcar2(onda, box, offset, S, sigma0)
wave onda
variable box //size of the boxcar
variable offset, S, sigma0

wavestats/Q onda
variable k=0, b=0, sd=0, i=0, l, m=0
variable length=(V_npnts)
make/O/N=(box) temp
variable tempb
variable tempb2

do	
	temp=onda[i+p]
	
	wavestats/Q temp
	k=k+ V_avg
	sd=sd+V_sdev
	//calculates the molecular brightness (epsilon) value for each box
	tempb=((V_sdev^2-sigma0^2)/(V_avg-offset))
	tempb2=b
	b=tempb2+tempb
		
i+=1
while(i<(length-box))
	
	sd=sd/i
	k=k/i
	b=b/i
	
make/O/N=3 out={k,sd,b}

wave wout=$("out")
return wout


end
/////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////Auxiliary Functions///////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////
/////////Polygon Form to select the brightness values within the polygon/////
/////////in the B vs I plot//////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////


function highlightBch1_poly()
wave intensity_cal=$("intensity_cal1")
wave brightness_cal=$("brightness_cal1")

nvar sizex, sizey
svar wn
variable getpoint, i,j

string ctrlName
variable V_npnts, V_stdev

make/O W_inpoly

make/O/N=(sizex*sizey) tempbrightness1
make/O/N=(sizex*sizey) brightness2Hist, intensity2avg
make/O/N=1 meanint, meanB

tempbrightness1=NaN

////this is to create the control panel

	DoWindow/F BfracS_vs_I_ch1						// Bring graph to front
	
	if (V_Flag == 0)									// Verify that graph exists
		Abort "UserCursorAdjust: No such graph."
		return -1
	endif

	NewPanel/K=2 /W=(139,341,382,432) as "Pause for Cursor"
	DoWindow/C tmp_PauseforCursor					// Set to an unlikely name
	AutoPositionWindow/E/M=1/R=BfracS_vs_I_ch1	// Put panel near the graph

	DrawText 21,20,"Adjust the polygon cursors and then. Doubleclick when finished"
	DrawText 21,40,"Click Continue."
	
	Button button1,pos={20,58},size={92,20}, title="Draw Polygon "
	Button button1, proc=puntinelpoligono

	Button button0,pos={120,58},size={92,20}, title="Continue"
	Button button0, proc=autochiuditi
	
	PauseForUser tmp_PauseforCursor, BfracS_vs_I_ch1
	
	Redimension /N=(sizex*sizey) intensity_cal, brightness_cal
	
	FindPointsInPoly intensity_cal, brightness_cal, x_poly, y_poly
	

variable k=0

for(i=0; i<(sizex*sizey); i+=1)
			
		if(W_inpoly[i]==1 && intensity_cal[i]!=NaN)//if a point has been selected by the cursors, use it for calculating the histogram of the brightness
			 tempbrightness1[i]=0
			 brightness2Hist[k]=brightness_cal[i]
			 intensity2avg[k]=intensity_cal[i]
			k+=1
		else
			tempbrightness1[i]=NaN
		endif
		
	
endfor

redimension/N=(sizex, sizey) tempbrightness1, intensity_cal, brightness_cal	

deletepoints (k), (sizex*sizey-k), brightness2Hist
deletepoints (k), (sizex*sizey-k), intensity2avg
meanint[0]=mean(intensity2avg)
meanB[0]=mean(brightness2Hist)
	
	wavestats brightness2Hist
	//do an histogram here
	brightness_Histogram("intensity2avg", "brightness2Hist")
//	

DoWindow/F brightness_ch1
RemoveImage tempbrightness1

appendimage/W=brightness_ch1 tempbrightness1
ModifyImage/W=brightness_ch1 tempbrightness1 ctab= {*,1,Grays,1}

print "mean Brightness ROI: "+num2str(meanB[0])
print "mean intensity ROI: "+num2str(meanint[0])

/////////
duplicate/O intensity2avg, $("I_"+wn)
duplicate/O brightness2hist, $("BR_"+wn)

////////
end

//////////////////////////////////////////////////////////////////////////////////////////////////
/////////Polygon form to select the brightness values within the polygon from intensity image/////
//////////////////////////////////////////////////////////////////////////////////////////////////

function highlightIch1_poly() //selects a region of interest in the intensity plot

wave intensity_cal=$("intensity_cal1")
wave brightness_cal=$("brightness_cal1")

nvar sizex, sizey
svar wn
variable getpoint, i,j

string ctrlName
variable V_npnts, V_stdev

make/O W_inpoly_I

duplicate/O intensity_cal, M_ROIMask
duplicate/O intensity_cal, intensity_roi
duplicate/O brightness_cal, brightness_roi
make/O/N=(sizex, sizey) tempbrightness1

//this is to create the control panel

	DoWindow/F intensity_ch1						// Bring graph to front
	
	if (V_Flag == 0)									// Verify that graph exists
		Abort "UserCursorAdjust: No such graph."
		return -1
	endif

	NewPanel/K=2 /W=(139,341,382,432) as "Pause for Cursor"
	DoWindow/C tmp_PauseforCursor_I					// Set to an unlikely name
	AutoPositionWindow/E/M=1/R=intensity_ch1		// Put panel near the graph

	DrawText 21,20,"Adjust the polygon cursors and then. Doubleclick when finished"
	DrawText 21,40,"Click Continue."
	
	Button button1,pos={20,58},size={92,20}, title="Draw Polygon "
	Button button1, proc=puntinelpoligono_I

	Button button0,pos={120,58},size={92,20}, title="Continue"
	Button button0, proc=autochiuditi_I
	
	PauseForUser tmp_PauseforCursor_I, intensity_ch1
	
//put here something to wait for user input.
	
	ImageBoundaryToMask width=sizex, height=sizey, xwave=x_poly, ywave=y_poly, seedX=mean(x_poly), seedY=mean(y_poly) 

	redimension/N=(sizex*sizey) intensity_roi, brightness_roi

variable k=0

	for(i=0; i<(sizex); i+=1)
		
		for (j=0; j<sizey; j+=1)
					
			if(M_ROIMask[i][j]==1)//if a point has been selected by the cursors, use it for calculating the histogram of the brightness
				 tempbrightness1[i][j]=1
				 intensity_roi[k]=intensity_cal[i][j]
				 brightness_roi[k]=brightness_cal[i][j]
				 k+=1
			else
				tempbrightness1[i][j]=NaN
			endif
			
		endfor		
	
	endfor

	
	deletepoints (k), (sizex*sizey-k), intensity_roi
	deletepoints (k), (sizex*sizey-k), brightness_roi
		
	brightness_Histogram("intensity_roi", "brightness_roi")
	
	DoWindow/F brightness_ch1
	appendimage/W=brightness_ch1 tempbrightness1
	ModifyImage/W=brightness_ch1 tempbrightness1 ctab= {*,1,Grays,1}
	DoWindow/F intensity_ch1
	RemoveFromGraph y_poly
	
	duplicate/O brightness_roi, $("BR_"+wn)
	duplicate/O intensity_roi, $("I_"+wn)
	
end


/////////////////////////////////////////////////////////////////////

function brightness_Histogram(intensity2avg_s, brightness2Hist_s)
string intensity2avg_s, brightness2Hist_s

wave intensity2avg_w=$(intensity2avg_s)
wave brightness2Hist_w=$(brightness2Hist_s)
variable NpuntiHist
variable V_stdev, V_npnts
	
	//do an histogram here
	wavestats/Q brightness2Hist_w
	NpuntiHist=3.49*V_stdev*(V_npnts)^(-1/3)
		
		if(NpuntiHist<=1)
		NpuntiHist=200
		Make/N=(NpuntiHist) /O W_Hist;DelayUpdate
		SetScale/P x 0,0.01,"", W_Hist
		else
		Make/N=(NpuntiHist) /O W_Hist;DelayUpdate
		endif
		
		print NpuntiHist
	
	Histogram/B=2 brightness2Hist_w,W_Hist;DelayUpdate
	
	dowindow/K B_Hist_ch1
	display/N=B_Hist_ch1 /M/W=(15,15,30,25 ); appendtograph/W=B_Hist_ch1 W_Hist
	ModifyGraph fSize=16,axThick=2
	Label left "Occurrences"
	Label bottom "Brightness"
	ShowInfo


				
				wavestats/Q intensity2avg_w //calculates the statistics of the pixels of the intensity image within the polygonal ROI
		    	make/O/N=(V_npnts) theory_data
		    	wave W_KSResults
			 	theory_data=gnoise(V_sdev/sqrt(2)) + V_avg //generates Gaussian distributed intensity values for a ROI containing as many pixels as pixelsinROI and its mean and variance
		 		StatsKSTest/ALPH=0.05/Z/Q/T=1 intensity2avg_w,theory_data //Calculates for Gaussianity using Komogorov-Smirnov criterion
	 	 		print "the gaussianity score is: "+num2str(W_KSResults[5]/W_KSResults[4]) //the bigger this value the better (it is the ratio between threshold and measured value, the measured value should be below threshold for gaussianity


print "mean Brightness ROI: "+num2str(V_avg)
print "mean intensity ROI: "+num2str(mean(intensity2avg_w))

end



////////////////////////////////
//////Button Functions//////////
////////////////////////////////


function autochiuditi(ctrlName): ButtonControl
String ctrlName

Dowindow/K tmp_PauseforCursor

end

/////////////////////////

function autochiuditi_I(ctrlName): ButtonControl
String ctrlName

Dowindow/K tmp_PauseforCursor_I

end

/////////////////////////

function autochiuditi_preview(ctrlName): ButtonControl
string ctrlName

Dowindow/K tmp_PauseforCursor_preview


end

/////////////////////////

function puntinelpoligono(ctrlName) : ButtonControl
String ctrlName

Dowindow/F BfracS_vs_I_ch1 

GraphWaveDraw/W=BfracS_vs_I_ch1 /O y_poly, x_poly

end

/////////////////////////

function puntinelpoligono_I(ctrlName) : ButtonControl
String ctrlName

Dowindow/F intensity_ch1 

GraphWaveDraw/W=intensity_ch1 /O y_poly, x_poly

end

/////////////////////////

function puntinelpoligono_preview(ctrlName) : ButtonControl
String ctrlName

Dowindow/F average_intensity_preview 

GraphWaveDraw/W=average_intensity_preview /O y_poly, x_poly

end

/////////////////////////

function autochiuditi_SpIDA_preview(ctrlName): ButtonControl
string ctrlName

Dowindow/K tmp_PauseforCursor_preview


end


////////////////////////

function puntinelpoligono_SpIDA_preview(ctrlName) : ButtonControl
String ctrlName

Dowindow/F SpIDA_preview

GraphWaveDraw/W=SpIDA_preview /O y_poly, x_poly

end
