#!/usr/bin/env wolframscript
(* ::Package:: *)

(* Load VMM Package *)
Needs["VMM`"]

(* Helper Functions *)
WriteLineToStdOut[str_]:= WriteString["stdout",str<>"\n"];

GetImageFiles[dir_]:= FileNames["*.png"|"*.jpg",{dir}];

ParseInputs[]:= Module[{args,inputScriptFile,inputDir,imgList,configFile,configAssoc},
	args = $ScriptCommandLine;
	inputDir = "";
	imgList = {};
	configFile = "";
	
	If[Length[args] == 3, 
		(* Correct number of args *)
		{inputScriptFile,inputDir,configFile}=args;
		
		If[DirectoryQ[inputDir],
			imgList=GetImageFiles[inputDir];,
			WriteLineToStdOut["Input argument " <> inputDir <> " could not be found"];
		];
		
		If[FileExistsQ[configFile],
			configAssoc=Import[configFile,"RawJSON"];,
			WriteLineToStdOut["Input argument " <> configFile <> " could not be found"];
		];
		,
		
		(* Incorrect number of args *)
		WriteLineToStdOut["Wrong number of input arguments recieved!"]
	];
	
	Return[{inputDir, imgList, configFile, configAssoc}];
];

(* Main *)
WriteLineToStdOut["\n\n**      ** ****     **** ****     ****
/**     /**/**/**   **/**/**/**   **/**
/**     /**/**//** ** /**/**//** ** /**
//**    ** /** //***  /**/** //***  /**
 //**  **  /**  //*   /**/**  //*   /**
  //****   /**   /    /**/**   /    /**
   //**    /**        /**/**        /**
    //     //         // //         // \n"]
WriteLineToStdOut["~~~~~~~~~~~~ VMM Software 0.3 ~~~~~~~~~~~"]
WriteLineToStdOut["Time to create an awesome microscopic mosaic!\n"]

(* Parse input arguments *)
{inputDir, imgList, configFile, configAssoc} = ParseInputs[];

If[Length[imgList]>50,
	WriteLineToStdOut["Input directory is: "<>inputDir];
	WriteLineToStdOut[("Processing " <> ToString[Length[imgList]] <> " image files found in " <> inputDir<>"\n")];
	proceed=True,
	WriteLineToStdOut["Too few images input! Please input at least 50 images."];
	proceed=False
];

WriteLineToStdOut["Configuration parameters found in "<> configFile];
WriteLineToStdOut[ToString[Column[Normal[configAssoc]]]<>"\n"];
configParams = {"downsample_factor","feature_dectector","export_intermediates","perform_inpaint","perform_brightness_equalization","qc_registered_shifts"};
If[ContainsAll[Keys[configAssoc],configParams],
	downsampleOk = IntegerQ[configAssoc[["downsample_factor"]]] && configAssoc[["downsample_factor"]]>=1;
	featureDetectorOk = ContainsAny[{configAssoc[["feature_dectector"]]},{"SURF","AKAZE","KAZE","ORB"}];
	boolInputsOk = And@@(BooleanQ/@Values[configAssoc[[{"export_intermediates","perform_inpaint","perform_brightness_equalization","qc_registered_shifts"}]]]);
	If[And@@{downsampleOk,featureDetectorOk,boolInputsOk},
		proceed=True,
		WriteLineToStdOut["One or more of the inputs from the configuration file was not successfully parsed!"];
		If[Not@downsampleOk,WriteLineToStdOut["downsample_factor should be a positive integer"]];
		If[Not@featureDetectorOk,WriteLineToStdOut["feature_dectector should be one of the following: " <> ToString[{"SURF","AKAZE","KAZE","ORB"}]]];
		If[Not@boolInputOk,
			boolInputsStr = StringJoin[Riffle[{"export_intermediates","perform_inpaint","perform_brightness_equalization","qc_registered_shifts"}," "]];
			WriteLineToStdOut[boolInputsStr<>" should all be Boolean inputs (True/False)"];
		];
		proceed=False
	];,
	WriteLineToStdOut["Config file missing one of the following parameters: " <> ToString[configParams]];
	proceed=False
];

If[proceed,
	(* Example *)
	(* VMM[inputDir, imgList, configFile, 1, "SURF", True, True, True, True] *)
	VMM[inputDir, 
		imgList, 
		configFile, 
		configAssoc[["downsample_factor"]],  
		configAssoc[["feature_dectector"]], 
		configAssoc[["export_intermediates"]], 
		configAssoc[["perform_inpaint"]], 
		configAssoc[["perform_brightness_equalization"]],
		configAssoc[["qc_registered_shifts"]]
	],
	WriteLineToStdOut["\nPlease call this program using:\n\n\t\t 'wolframscript -file vmm.wls DIR CONFIG_FILE' \n\nand pass two arguments: 1) a directory containing image files (.jpg or .png) and 2) a config.json with other parameters \n"]
];






