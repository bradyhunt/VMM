(* ::Package:: *)

BeginPackage["VMM`"]

VMM::usage =
		"VMM[inputDir, imgFileList, configFile, exludeIntermediates, downsampleFactor, featureDetector, inpaint, brigthnessEqualization, applyQC] generates and writes out a mosaicked image using the input image files."

Begin["`Private`"]


MeanImage[imgs_] := Module[{imgType,maxValue,total,normalize,mean},
	imgType = ImageType[First@imgs];
	total = Total[Map[(ImageData[#,imgType]&),imgs]];
	normalize = ConstantArray[1/(Length[imgs]),Dimensions[First[imgs]]];
	mean = Image[total*normalize,ImageType[First[imgs]]];	
	Return[mean];
];

(* ROI Mask Functions *)
EdgeCoords[binaryMask_]:=Position[Reverse/@Transpose[IntegerPart[ImageData[binaryMask]]],1];
FitCircleToPoints[pts_]:=Module[{reg,lm,bf,exp,center,rad,fitError,fitAssoc},
	reg={2 #1,2 #2,#2^2+#1^2}&@@@pts;
	lm=LinearModelFit[reg,{1,x,y},{x,y}];
	bf=lm["BestFitParameters"];
	exp=(x-#2)^2+(y-#3)^2-#1-#2^2-#3^2&@@bf;
	{center,rad}={{#2,#3},Sqrt[#2^2+#3^2+#1]}&@@bf;
	fitError=Mean[Abs[lm["FitResiduals"]]];
	fitAssoc=Association[{"expression"->exp,"center"->center,"radius"->rad,"error"->fitError}];
	Return[fitAssoc];
];

CircFit[edgeCoords_,numPoints_]:=Module[{edgeCoordsReduced,fit},
	SeedRandom[7];
	edgeCoordsReduced = If[Length[edgeCoords]>numPoints,RandomSample[edgeCoords,UpTo[numPoints]],edgeCoords];
	fit = FitCircleToPoints[edgeCoordsReduced];
	Return[fit];
];

CircleROIMask[img_,center_,radius_]:= Module[{bg,disk,mask},
	bg=Image[ConstantArray[0,Reverse[ImageDimensions[img]]]];
	disk = Disk[center,radius];
	mask = Erosion[ImageResize[Binarize[Rasterize[Show[bg,Graphics[{White,disk}],ImageSize->ImageDimensions[img]/2]]],ImageDimensions[img]],1];
	Return[mask];
];

BinarizeMeanImage[meanImage_,blurRadius_:10]:= SelectComponents[Binarize[BrightnessEqualize[meanImage,GaussianFilter[meanImage,blurRadius]],Method->"Cluster"],"Area",-1];

MeanImageROIMask[meanImage_,innerRadiusCrop_:15]:=Module[{binaryMask,edgeMask,edgeCoords, roiMask,circfit,roiMaskInner,deadPixelMask},
	binaryMask = BinarizeMeanImage[meanImage];
	edgeMask = Dilation[EdgeDetect[FillingTransform[binaryMask]],1];
	circfit = CircFit[EdgeCoords[edgeMask],100];
	roiMask = CircleROIMask[meanImage,circfit["center"],circfit["radius"]];
	roiMaskInner = CircleROIMask[meanImage,circfit["center"],circfit["radius"]-innerRadiusCrop];
	deadPixelMask = ColorNegate[binaryMask]-Dilation[ColorNegate[roiMask],DiskMatrix[3]];
	Return[{binaryMask,roiMask,roiMaskInner,deadPixelMask}];
];

VMMBrightnessEqualize[img_,imgROIMask_,meanImg_]:= ImageMultiply[Lighter[Image[BrightnessEqualize[img,meanImg,ColorNegate[imgROIMask]],"Byte"]],imgROIMask];

WriteLineToStdOut[str_]:= WriteString["stdout",str<>"\n"];
FormatTimeString[t_]:=ToString[NumberForm[t,{4,3}]];

ExportFile[fileName_,expr_, outputDir_, exportFlag_]:= If[exportFlag,Export[FileNameJoin[{outputDir,fileName}],expr],OverwriteTarget->True];
ExportFiles[fileNames_,exprs_, outputDir_, exportFlag_]:=MapThread[(ExportFile[#1,#2,outputDir,exportFlag]&),{fileNames,exprs}];
ExportConfig[configFileInput_, outputDir_, exportFlag_]:= If[exportFlag, CopyFile[configFileInput, FileNameJoin[{outputDir,FileNameTake[configFileInput]}],OverwriteTarget->True]];

ImageSharpness[img_,imgROIMask_:Missing[],kernelSize_:2]:= Module[{imgSatMask,imgROIMask2,imgROIMask3,sharpnessMetric},
imgROIMask2 = If[MissingQ[imgROIMask],ConstantArray[1,Reverse[ImageDimensions[img]]],imgROIMask];
imgSatMask = Dilation[Binarize[img,0.9],DiskMatrix[5]];
imgROIMask3 = imgROIMask2-imgSatMask;
sharpnessMetric = ImageMeasurements[LaplacianGaussianFilter[img,kernelSize],"StandardDeviation",Masking->imgROIMask3];
Return[sharpnessMetric];
];

ImageShiftWithROI[imgs_,roiMask_,features_]:=Module[{img1,img2,matches,numPoints,err,tr,shiftX,shiftY,shiftNorm,returnAssoc},
	{img1,img2} = imgs;
	matches = ImageCorrespondingPoints[img1,img2,TransformationClass->"Translation",Method->features,Masking->roiMask];
	numPoints = Length[Transpose[matches]];
	If[numPoints==0, Return[AssociationThread[{"NUM_PTS","ERR","SHIFT_X","SHIFT_Y","SHIFT_NORM"}->{numPoints,Missing[],Missing[],Missing[],Missing[]}]]];
	{err,tr}=(FindGeometricTransform[#1,#2,TransformationClass->"Translation"]&)@@matches;
	{shiftX,shiftY}=Round[tr@{0,0}];
	shiftNorm = N[Norm[{shiftX,shiftY}]];
	returnAssoc=AssociationThread[{"NUM_PTS","ERR","SHIFT_X","SHIFT_Y","SHIFT_NORM"}->{numPoints,err,shiftX,shiftY,shiftNorm}];
	Return[returnAssoc];
];

RegisterShifts[images_,roiMask_,features_:"SURF"]:= Module[{shiftAssocs},
	
	shiftAssocs = ParallelMap[(ImageShiftWithROI[#,roiMask,features]&),Partition[images,2,1]];
	shiftAssocs = Prepend[shiftAssocs,<|"NUM_PTS"->Missing[],"ERR"->Missing[],"SHIFT_X"->0,"SHIFT_Y"->0,"SHIFT_NORM"->0.`|>];
	shiftAssocs[[All,"FRAME_NUM"]]=Range[Length[shiftAssocs]];
	shiftAssocs = shiftAssocs[[All,{"FRAME_NUM","NUM_PTS","ERR","SHIFT_X","SHIFT_Y","SHIFT_NORM"}]];
	
	Return[shiftAssocs];
];

QCShifts[images_, roiMaskInner_, shiftAssocs_]:=Module[{shiftAssocsCopy, neighborDistances,shiftConversionFactor,moveAvgSampleSize,
	centerCoord,registerErrorFunc,distanceParams,distanceParamsRandShuffled,shiftErrors,shiftErrorsRand},
	
	shiftAssocsCopy = shiftAssocs;
	neighborDistances=Append[Prepend[(N[Mean[{EuclideanDistance[#1,#2],EuclideanDistance[#2,#3]}]]&)@@@Partition[Values[shiftAssocsCopy[[All,{"SHIFT_X","SHIFT_Y"}]]]/.{Missing[]->0},3,1],0],0];
	shiftAssocsCopy[[All,"NEIGHBOR_DIST"]] = neighborDistances;
	shiftAssocsCopy[[All,"SHARPNESS"]]=ParallelMap[(ImageSharpness[#,roiMaskInner]&),images];
	
	distanceParams = MapThread[Append,{Partition[images,2,1],Rest[Values[shiftAssocsCopy[[All,{"SHIFT_X","SHIFT_Y"}]]]]/.{Missing[]->0}}];
	SeedRandom[7];
	distanceParamsRandShuffled = MapThread[Append,{Partition[RandomSample[images],2,1],Rest[Values[shiftAssocsCopy[[All,{"SHIFT_X","SHIFT_Y"}]]]]/.{Missing[]->0}}];
	
	centerCoord = Round[ImageDimensions[First@images]/2];
	SetSharedVariable[centerCoord];
	registerErrorFunc = (ImageDistance[#[[1]],#[[2]],centerCoord,If[Or@@MapThread[Greater,{Abs/@#[[3]],centerCoord}],centerCoord,centerCoord-#[[3]]],Masking->{roiMaskInner,roiMaskInner},DistanceFunction->"RootMeanSquare"]&);
	shiftErrors = ParallelMap[registerErrorFunc,distanceParams];
	shiftErrorsRand = ParallelMap[registerErrorFunc,distanceParamsRandShuffled];

	shiftAssocsCopy[[All,"SHIFT_RMSE"]] = Prepend[shiftErrors,Mean[shiftErrors]];
	shiftAssocsCopy[[All,"SHIFT_RMSE_SHUFFLE"]]=Prepend[shiftErrorsRand,Mean[shiftErrorsRand]];
	shiftAssocsCopy[[All,"SHIFT_RMSE_QC_10_PERCENT"]] = Quantile[shiftErrorsRand,0.1];
	shiftAssocsCopy[[All,"SHIFT_RMSE_QC_01_PERCENT"]] = Quantile[shiftErrorsRand,0.01];
	shiftAssocsCopy[[All,"SHIFT_RMSE_QC"]] =(If[#"SHIFT_RMSE">#"SHIFT_RMSE_QC_01_PERCENT"||#"SHARPNESS"<0.009,1,0]& )/@shiftAssocsCopy;
	
	Return[shiftAssocsCopy];
];

SetAttributes[UpdateMosaicInPlaceMask,HoldFirst];
UpdateMosaicInPlaceMask[mosaic_,topLeftCoord_,newFrame_,newFramSubtract_,oldFrameSubstract_]:=Module[{frameColSize,frameRowSize,col,row,oldFrame,newFrameMasked,oldFrameMasked,frameUpdate},
	{frameColSize,frameRowSize} = Dimensions[newFrame];
	{col,row} = topLeftCoord;
	oldFrame=mosaic[[col;;(col+frameColSize-1),row;;(row+frameRowSize-1)]];
	newFrameMasked = newFrame-newFramSubtract;
	oldFrameMasked = oldFrame-oldFrameSubstract;
	frameUpdate=Map[Max,(Join[##,3]&)@@((ArrayReshape[#,Append[Dimensions[#],1]]&)/@{oldFrameMasked,newFrameMasked}),{2}];
	mosaic[[col;;(col+frameColSize-1),row;;(row+frameRowSize-1)]]=frameUpdate;
];

BuildMosaic[images_,shiftAssocs_,roiMaskInner_, excludeQCFrames_] := Module[
	{shiftAssocsQC,padX,padY,pathNoQC,path,plotRange,pathPlot,mosaicDims,w,h,
	mosaic,mosaicPathCoords,imgType,framesData,frameNewSubtract,frameOldSubtract},
	
	padX=250;
	padY=250;
	
	shiftAssocsQC = shiftAssocs;
	If[excludeQCFrames, shiftAssocsQC[[Flatten[Position[shiftAssocs[[All,"SHIFT_RMSE_QC"]],1]],{"SHIFT_X","SHIFT_Y","SHIFT_NORM"}]]={Missing[],Missing[],Missing[]}];
	path = Accumulate[Values[shiftAssocsQC[[All,{"SHIFT_X","SHIFT_Y"}]]]/.{Missing[]->0}];
	path = (#+(Abs/@Min/@Transpose[path])&)/@path;
	plotRange = Transpose[{Floor[(Min/@Transpose[path])-padX,100],Ceiling[(Max/@Transpose[path])+padY,100]}];
	pathPlot = Rasterize[ListPlot[path,PlotRange->plotRange,AspectRatio->Automatic],ImageResolution->144];
	
	mosaicDims = Abs/@Subtract@@@plotRange;
	{w,h} = mosaicDims;
	mosaic = ConstantArray[0,Reverse@mosaicDims];
	mosaicPathCoords = (({-padY,padX})+{-1,1}+#&)/@({h-#1,#2}&)@@@ Reverse/@ path;
	imgType = ImageType[First@images];
	framesData = (ImageData[#,imgType]&)/@images;
	frameNewSubtract=ImageData[ColorNegate[roiMaskInner],imgType];
	frameOldSubtract=ImageData[roiMaskInner,imgType];
	
	Do[UpdateMosaicInPlaceMask[mosaic,mosaicPathCoords[[n]],framesData[[n]],frameNewSubtract,frameOldSubtract],{n,1,Length[framesData]}];
	Return[{pathPlot,Image[mosaic,imgType]}];
];

VMM[inputDir_, imgList_, configFile_, downsample_:1, features_:"SURF", exportIntermediates_:True, performInpaint_:True, performBrightnessEqualization_:True, qcRegistrations_:True]  := Module[
	{imgListCopy, imgListProcessed, processedDir, kernels,outputDir,images,meanImage,masks,binaryMask,roiMask,roiMaskInner,deadPixelMask,maskFiles,
	shiftAssocs,pathPlot,mosaic,mosaicFile,exportPath,blank,t0,t1,t2,t3,t4,t5,t6,t7,t8,t9,timeTotal},
	
	If[downsample!=1,WriteLineToStdOut["Downsampling input frames by factor of " <> ToString[downsample]]];
	imgListCopy = Downsample[imgList,downsample];
	outputDir = inputDir<>"_vmm_mosaic";
	If[Not[DirectoryQ[outputDir]],CreateDirectory[outputDir]];
	WriteLineToStdOut["Exporting output files to: " <> outputDir];
	ExportConfig[configFile,outputDir,True];
	
	WriteString["stdout", "VMM [0/9] - Launching kernels... "];
	{t0,kernels} = LaunchKernels[] // AbsoluteTiming;
	WriteString["stdout", "\t\t("<>FormatTimeString[t0]<>" seconds)\n"];
	
	WriteString["stdout", "VMM [1/9] - Importing images... "];
	{t1, images} = ParallelMap[Import,imgListCopy] // AbsoluteTiming;
	WriteString["stdout", "\t\t("<>FormatTimeString[t1]<>" seconds)\n"];
	
	WriteString["stdout", "VMM [2/9] - Computing mean image... "];
	{t2,meanImage} = MeanImage[images] // AbsoluteTiming;
	ExportFile["mean_image.png", meanImage,outputDir,exportIntermediates];
	WriteString["stdout", "\t\t("<>FormatTimeString[t2]<>" seconds)\n"];
	
	WriteString["stdout", "VMM [3/9] - Creating ROI mask... "];
	{t3,masks} = MeanImageROIMask[meanImage] // AbsoluteTiming;
	{binaryMask,roiMask,roiMaskInner,deadPixelMask} = masks;
	maskFiles = {"binary_mask.png","roi_mask.png","roi_inner_mask.png","dead_pixel_mask.png"};
	ExportFiles[maskFiles,masks,outputDir,exportIntermediates];
	WriteString["stdout", "\t\t("<>FormatTimeString[t3]<>" seconds)\n"];
	
	WriteString["stdout", "VMM [4/9] - Interpolate static pixels... "];
	If[performInpaint,
	SetSharedVariable[deadPixelMask];
	{t4,images} = ParallelMap[(Inpaint[#,deadPixelMask]&),images] // AbsoluteTiming;
	meanImage = Inpaint[meanImage,deadPixelMask];
	ExportFile["mean_image_2.png",meanImage,outputDir,exportIntermediates];
	WriteString["stdout", "\t("<>FormatTimeString[t4]<>" seconds)\n"];,
	t4=0;WriteString["stdout", "\t SKIPPED \n"];
	];
	
	WriteString["stdout", "VMM [5/9] - Brightness equalization... "];
	If[performBrightnessEqualization,
	SetSharedVariable[meanImage,roiMask];
	{t5,images} = ParallelMap[(VMMBrightnessEqualize[#,roiMask,meanImage]&),images] // AbsoluteTiming;
	WriteString["stdout", "\t\t("<>FormatTimeString[t5]<>" seconds)\n"];,
	t5=0;WriteString["stdout", "\t\t SKIPPED \n"];
	];
	
	WriteString["stdout", "VMM [6/9] - Export processed frames... "];
	If[exportIntermediates,
	{t6,blank}={processedDir = FileNameJoin[{outputDir,FileBaseName[inputDir]<>"_processed"}];
	imgListProcessed = (FileNameJoin[{FileBaseName[inputDir]<>"_processed",#}]&)/@FileNameTake/@imgListCopy;
	If[exportIntermediates, If[Not[DirectoryQ[processedDir]],CreateDirectory[processedDir]]];
	ExportFiles[imgListProcessed,images,outputDir,exportIntermediates];}//AbsoluteTiming;
	WriteString["stdout", "\t\t("<>FormatTimeString[t6]<>" seconds)\n"],
	t6=0;WriteString["stdout", "\t\t SKIPPED \n"];
	];
	
	WriteString["stdout", "VMM [7/9] - Frame registration... "];
	{t7,shiftAssocs} = RegisterShifts[images,roiMask,features] // AbsoluteTiming;
	WriteString["stdout", "\t\t("<>FormatTimeString[t7]<>" seconds)\n"];
	
	WriteString["stdout", "VMM [8/9] - QC registrations... "];
	If[qcRegistrations,
	{t8,shiftAssocs} = QCShifts[images,roiMaskInner,shiftAssocs] // AbsoluteTiming;
	WriteString["stdout", "\t\t("<>FormatTimeString[t8]<>" seconds)\n"];,
	t8=0;WriteString["stdout", "\t\t SKIPPED \n"];
	];
	
	WriteString["stdout", "VMM [9/9] - Building mosaic... "];
	{t9,{pathPlot,mosaic}} = BuildMosaic[images,shiftAssocs,roiMaskInner,qcRegistrations] // AbsoluteTiming;
	(* always export registrations, path plot, and final mosaic *)
	ExportFile["vmm_shifts.csv",Dataset@shiftAssocs,outputDir,True];
	ExportFile["path_plot.png",pathPlot,outputDir,True];
	ExportFile["mosaic.png", mosaic, outputDir, True];
	WriteString["stdout", "\t\t\t("<>FormatTimeString[t9]<>" seconds)\n"];
	
	WriteString["stdout", "\nGenerated mosaic saved to "<>FileNameJoin[{outputDir,"mosaic.png"}]<>"!\n"];
	timeTotal = Total[{t0,t1,t2,t3,t4,t5,t6,t7,t8,t9}];
	WriteString["stdout", "Total time: "<>FormatTimeString[timeTotal]<>"\n"];
];

End[ ]

EndPackage[ ]

