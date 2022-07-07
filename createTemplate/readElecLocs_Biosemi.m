function [fifIdx, elecIdx, elecCoord] = readElecLocs(eegSystem)

fifIdx = [ 22 9 11 33 24 36 45 52 58 62 70 75 83 92 96 108 104 124 122];
% matching electrodes between systems
matchElec = {'Fp1' 'Fp2' 'Fz' 'F7' 'F3' 'C3' 'T7' 'P3' 'P7' 'Pz' 'O1' 'Oz' 'O2' 'P4' 'P8' 'T8' 'C4' 'F4' 'F8'};
        
switch eegSystem
     
    case 'Standard_1005' 
        addpath (genpath('/Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork/'))
        A = ft_read_sens('EEG systems/standard/standard_1005.elc');
        elecCoord = A.elecpos;
        elecIdx = zeros(length(matchElec),1);
        for mm = 1:length(matchElec)
            elecIdx(mm) = find(strcmp(A.label,matchElec(mm)));
        end
        
    case 'EGI32' 
        elecIdx= [1 2 17 11 3 5 13 7 17 19 9 20 10 8 16 14 6 4 12];        
        
        fid = fopen('EEG systems/EGI/AdultAverageNet32_v1.sfp');
        A = textscan(fid,'%*s %f %f %f'); % skip 1st column of names
        A = cell2mat(A);
        elecCoord = A(4:end-1,:); % skip 1st 3 fiducials & Cz
        fclose(fid);  

    case 'EGI64' 
        elecIdx= [ 10 5 6 18 12 20 24 28 30 34 35 37 39 42 44 52 50 60 58];
        
        fid = fopen('EEG systems/EGI/AdultAverageNet64_v1.sfp');
        A = textscan(fid,'%*s %f %f %f'); % skip 1st column of names
        A = cell2mat(A);
        elecCoord = A(4:end-1,:); % skip 1st 3 fiducials & Cz
        fclose(fid);
        
    case 'EGI128' 
        elecIdx= fifIdx;
%         fid = fopen('EEG systems/EGI/GSN-HydroCel-129.sfp');
        fid = fopen('EEG systems/EGI/AdultAverageNet128_v1.sfp');
        A = textscan(fid,'%*s %f %f %f'); % skip 1st column of names
        A = cell2mat(A);
        elecCoord = A(4:end-1,:); % skip 1st 3 fiducials & Cz
        fclose(fid);
    
    case 'EGI256'
        elecIdx= [37 18 21 47 36 59 69 87 96 101 116 126 150 153 170 202 183 224 2];
        
        fid = fopen('EEG systems/EGI/AdultAverageNet256_v1.sfp');
        A = textscan(fid,'%*s %f %f %f'); % skip 1st column of names
        A = cell2mat(A);
        elecCoord = A(4:end-1,:); % skip 1st 3 fiducials & Cz
        fclose(fid);
        
    case 'biosemi32' % using circumference = 55 cm
        elecIdx= [ 1 30 31 3 4 8 7 12 11 13 15 16 17 19 20 24 23 27 28];
        elecCoord = .001*[ -27	83	-3
            -36	76	24
            -71	51	-3
            -48	59	44
            -33	33	74
            -78	30	27
            -87	0	-3
            -63	0	61
            -33	-33	74
            -78	-30	27
            -71	-51	-3
            -48	-59	44
            0	-63	61
            -36	-76	24
            -27	-83	-3
            0	-87	-3
            27	-83	-3
            36	-76	24
            48	-59	44
            71	-51	-3
            78	-30	27
            33	-33	74
            63	0	61
            87	0	-3
            78	30	27
            33	33	74
            48	59	44
            71	51	-3
            36	76	24
            27	83	-3
            0	63	61
            0	0	88];
        
    case 'biosemi64'      
        elecIdx= [ 1 34 38 7 5 13 15 21 23 31 27 29 64 58 60 52 50 40 42];
        
        elecCoord = .001*[ -27	83	-3
            -51	71	-3
            -36	76	24
            -25	62	56
            -48	59	44
            -64	55	23
            -71	51	-3
            -83	27	-3
            -78	30	27
            -59	31	56
            -33	33	74
            -34	0	81
            -63	0	61
            -82	0	31
            -87	0	-3
            -83	-27	-3
            -78	-30	27
            -59	-31	56
            -33	-33	74
            -25	-62	56
            -48	-59	44
            -64	-55	23
            -71	-51	-3
            -64	-47	-37
            -51	-71	-3
            -36	-76	24
            -27	-83	-3
            0	-79	-37
            0	-87	-3
            0	-82	31
            0	-63	61
            0	-34	81
            0	87	-3
            27	83	-3
            51	71	-3
            36	76	24
            0	82	31
            0	63	61
            25	62	56
            48	59	44
            64	55	23
            71	51	-3
            83	27	-3
            78	30	27
            59	31	56
            33	33	74
            0	34	81
            0	0	88
            34	0	81
            63	0	61
            82	0	31
            87	0	-3
            83	-27	-3
            78	-30	27
            59	-31	56
            33	-33	74
            25	-62	56
            48	-59	44
            64	-55	23
            71	-51	-3
            64	-47	-37
            51	-71	-3
            36	-76	24
            27	-83	-3];
        
    case 'biosemi128'
        elecIdx= [ 93    80    85   103   100   115   119     7   127    19    15 ...
            23    28    36    43    58    54    68    71];
        
        elecCoord = .001*[ 0	0	88
            0	-17	86
            0	-34	81
            0	-50	72
            -24	-58	61
            -45	-45	61
            -52	-52	47
            -48	-66	31
            -51	-70	14
            -51	-71	-3
            -50	-69	-20
            -47	-64	-37
            -25	-75	-37
            -26	-81	-20
            -27	-83	-3
            -27	-82	14
            -25	-78	31
            -28	-68	47
            0	-63	61
            0	-74	47
            0	-82	31
            0	-86	14
            0	-87	-3
            0	-85	-20
            0	-79	-37
            25	-75	-37
            26	-81	-20
            27	-83	-3
            27	-82	14
            25	-78	31
            28	-68	47
            24	-58	61
            17	-5	86
            24	-24	81
            45	-45	61
            52	-52	47
            48	-66	31
            51	-70	14
            51	-71	-3
            50	-69	-20
            47	-64	-37
            69	-50	-20
            71	-51	-3
            70	-51	14
            66	-48	31
            83	-27	-3
            82	-27	14
            78	-25	31
            68	-28	47
            58	-24	61
            43	-25	72
            34	0	81
            50	0	72
            63	0	61
            74	0	47
            82	0	31
            86	0	14
            87	0	-3
            83	27	-3
            82	27	14
            78	25	31
            68	28	47
            58	24	61
            43	25	72
            10	14	86
            24	24	81
            45	45	61
            52	52	47
            66	48	31
            70	51	14
            71	51	-3
            51	71	-3
            51	70	14
            48	66	31
            25	43	72
            24	58	61
            28	68	47
            25	78	31
            27	82	14
            27	83	-3
            0	87	-3
            0	86	14
            0	82	31
            0	74	47
            0	63	61
            0	50	72
            0	34	81
            -25	43	72
            -24	58	61
            -28	68	47
            -25	78	31
            -27	82	14
            -27	83	-3
            -51	71	-3
            -51	70	14
            -48	66	31
            -10	14	86
            -24	24	81
            -45	45	61
            -52	52	47
            -66	48	31
            -70	51	14
            -71	51	-3
            -83	27	-3
            -82	27	14
            -78	25	31
            -68	28	47
            -58	24	61
            -43	25	72
            -34	0	81
            -17	-5	86
            -24	-24	81
            -43	-25	72
            -50	0	72
            -63	0	61
            -74	0	47
            -82	0	31
            -86	0	14
            -87	0	-3
            -83	-27	-3
            -82	-27	14
            -78	-25	31
            -68	-28	47
            -58	-24	61
            -66	-48	31
            -70	-51	14
            -71	-51	-3
            -69	-50	-20];
        
    case 'biosemi256'
        elecIdx= [ 157 128 145 188 208 203 222 6 253 19 41 75 82 87 108];
        
        elecCoord = .001*[ 0	0	88
            0	-14	86
            0	-28	83
            0	-41	78
            0	-52	70
            0	-63	61
            -16	-61	61
            -19	-69	50
            -16	-77	38
            -17	-82	25
            -17	-85	11
            -17	-86	-3
            -17	-84	-17
            -17	-80	-31
            -16	-74	-43
            0	-76	-43
            0	-82	-31
            0	-86	-17
            0	-87	-3
            0	-87	11
            0	-84	25
            0	-79	38
            0	-72	50
            16	-61	61
            19	-69	50
            16	-77	38
            17	-82	25
            17	-85	11
            17	-86	-3
            17	-84	-17
            17	-80	-31
            16	-74	-43
            16	-22	83
            16	-37	78
            31	-42	70
            31	-55	61
            36	-62	50
            32	-72	38
            34	-77	25
            33	-80	11
            33	-81	-3
            33	-79	-17
            33	-75	-31
            31	-70	-43
            48	-66	-31
            48	-71	-17
            49	-73	-3
            48	-72	11
            49	-68	25
            46	-64	38
            51	-51	50
            45	-45	61
            30	-27	78
            42	-31	70
            55	-31	61
            62	-36	50
            59	-53	38
            62	-56	25
            61	-61	11
            62	-62	-3
            61	-61	-17
            61	-55	-31
            71	-41	-31
            71	-48	-17
            13	-4	86
            26	-9	83
            39	-13	78
            50	-16	70
            61	-16	61
            69	-19	50
            75	-24	38
            68	-39	38
            73	-42	25
            72	-48	11
            73	-49	-3
            79	-33	-17
            81	-33	-3
            80	-33	11
            80	-26	25
            85	-17	11
            86	-17	-3
            87	0	-3
            87	0	11
            84	-9	25
            79	-8	38
            72	0	50
            63	0	61
            52	0	70
            40	4	78
            50	16	70
            61	16	61
            69	19	50
            79	8	38
            84	9	25
            85	17	11
            86	17	-3
            8	11	86
            26	9	83
            35	20	78
            42	31	70
            55	31	61
            75	24	38
            80	26	25
            80	33	11
            81	33	-3
            79	33	-17
            71	48	-17
            73	49	-3
            72	48	11
            73	42	25
            68	39	38
            62	36	50
            45	45	61
            51	51	50
            59	53	38
            62	56	25
            61	61	11
            62	62	-3
            61	61	-17
            49	73	-3
            48	72	11
            49	68	25
            46	64	38
            36	62	50
            32	72	38
            34	77	25
            33	80	11
            33	81	-3
            0	28	83
            16	22	83
            24	33	78
            31	42	70
            31	55	61
            16	61	61
            19	69	50
            16	77	38
            17	82	25
            17	85	11
            17	86	-3
            0	87	-3
            0	87	11
            0	84	25
            0	79	38
            0	72	50
            0	63	61
            0	52	70
            16	50	70
            8	40	78
            -8	40	78
            -16	50	70
            -16	61	61
            -19	69	50
            -16	77	38
            -17	82	25
            -17	85	11
            -17	86	-3
            -33	81	-3
            -33	80	11
            -34	77	25
            -32	72	38
            -8	11	86
            -16	22	83
            -24	33	78
            -31	42	70
            -31	55	61
            -36	62	50
            -46	64	38
            -49	68	25
            -48	72	11
            -49	73	-3
            -61	61	-17
            -62	62	-3
            -61	61	11
            -62	56	25
            -59	53	38
            -51	51	50
            -45	45	61
            -42	31	70
            -35	20	78
            -26	9	83
            -40	4	78
            -50	16	70
            -55	31	61
            -62	36	50
            -68	39	38
            -73	42	25
            -72	48	11
            -73	49	-3
            -71	48	-17
            -79	33	-17
            -81	33	-3
            -80	33	11
            -13	-4	86
            -26	-9	83
            -39	-13	78
            -52	0	70
            -61	16	61
            -69	19	50
            -75	24	38
            -80	26	25
            -85	17	11
            -86	17	-3
            -87	0	-3
            -87	0	11
            -84	9	25
            -79	8	38
            -72	0	50
            -63	0	61
            -69	-19	50
            -79	-8	38
            -84	-9	25
            -85	-17	11
            -86	-17	-3
            -79	-33	-17
            -81	-33	-3
            -80	-33	11
            -80	-26	25
            -75	-24	38
            -68	-39	38
            -73	-42	25
            -72	-48	11
            -73	-49	-3
            -71	-48	-17
            -71	-41	-31
            -16	-22	83
            -30	-27	78
            -42	-31	70
            -50	-16	70
            -61	-16	61
            -55	-31	61
            -62	-36	50
            -59	-53	38
            -62	-56	25
            -61	-61	11
            -62	-62	-3
            -61	-61	-17
            -61	-55	-31
            -48	-66	-31
            -48	-71	-17
            -49	-73	-3
            -48	-72	11
            -49	-68	25
            -46	-64	38
            -51	-51	50
            -45	-45	61
            -31	-42	70
            -16	-37	78
            -31	-55	61
            -36	-62	50
            -32	-72	38
            -34	-77	25
            -33	-80	11
            -33	-81	-3
            -33	-79	-17
            -33	-75	-31
            -31	-70	-43];
                
    otherwise
        disp('Sorry montage not currently available!')
end
