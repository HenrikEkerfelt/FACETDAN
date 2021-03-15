classdef DataSetDAN < handle
    %DATASETDAN Summary of this class goes here
    %   Detailed explanation goes here
    %
    %   Puclic properties:
    %       dataSet
    %       dataSetID
    %       hdr
    %
    %   Private properties:
    %       EPICS_INDEX   Indices of usable data
    %       epics_idx     Indices of currently used data
    %
    %       CAM_INDEX     Indices of usable data
    %       cam_idx       Indices of currently used data
    %
    %       isSorted
    %       sortLab
    %   
    %       visImageIncrement
    %
    %   Public methods:
    %       waterfallPlot
    %       correlationPlot
    %       visImages
    %
    %
    %   Private methods:
    %   
    
    properties (Access = public)
        dataSet;
        dataSetID;
        hdr;
    end
    
    properties (Access = private)
        EPICS_INDEX;
        CAM_INDEX;
        
        cam_idx;
        epics_idx;
        
        cam_select; 
        epics_select;
        
        cam_sort;
        epics_sort;
        
        include_data;
        
        isSorted = 0;
        sortLab = '';
        
        visImageIncrement = 1;
        
    end
    
    % Constructor
    methods 
        function s = DataSetDAN(dSID)
            %DATASETDAN Construct an instance of this class
            %   Detailed explanation goes here
            
            s.dataSetID = dSID;
            [s.dataSet,s.hdr] = getDataSet(dSID);
            disp('dataSet succefully loaded');
            s.select_cams();
            
            
        end
    end
        
    methods (Access = public) 
        
        function sortOnFSArray(s, FSArray)
        %% sortOnFSArray(FSArray) sorts the data for single diagnostic
        %   visualization on the scalar values in FSArray. 
        %   
        %   Example use:
        %
        %   dataSetDAN.sortOnFSArray({'BPMS_LI20_3315_X'})
        %       data is sorted on the scalar diagnostic 'BPMS_LI20_3315_X'
        %
        %   dataSetDAN.sortOnFSArray()
        %       resets to order based on shot number
        
            if nargin == 2
                type = s.hlpIsFSArray(FSArray);
                [FSArray,sortLab] = s.hlpGetScalarArray(FSArray, type);
                [~, I] = sort(FSArray);

                s.cam_sort = s.cam_idx(I);
                s.epics_sort = s.epics_idx(I);
                [~,idx] = intersect(s.cam_sort, s.cam_select);
                s.cam_idx = s.cam_sort(sort(idx));
                
                [~,idx2] = intersect(s.epics_select, s.epics_sort);
                s.epics_idx = s.epics_sort(sort(idx2));
                
                s.isSorted = 1;
                s.sortLab = sprintf('Sorted on: %s',sortLab)
                
            elseif nargin == 1
                s.cam_sort = s.CAM_INDEX;
                s.epics_sort = s.EPICS_INDEX;
                
                [~,idx] = intersect(s.cam_sort, s.cam_select);
                s.cam_idx = s.cam_sort(sort(idx));
                
                [~,idx2] = intersect(s.epics_sort, s.epics_select);
                s.epics_idx = s.epics_sort(sort(idx2));
                
                s.isSorted = 0;
                s.sortLab = '';
                
                
            end
            
        end
        
        function selectImages(s, FSArray, boolFcn)
        %% selectImages(s, FSArray, validFcn)
        % Of the already selected images, excludes the shots not fulfilling
        % the boolFcn.
        %
        % Example use:
        %
        % selectImages({'AX_IMG1', @(x) sum(sum(x)) }, @(x) x > 3e6)
        % selects the images from AX_IMG1 with a total pixel count sum 
        % greater than 3e6.
        % 
        % selectImages()
        % resets (removes) all previous selections.
        %
        
            if nargin == 3
                type = s.hlpIsFSArray(FSArray);
                [FSArray,FSLabel] = s.hlpGetScalarArray(FSArray, type);
                
                s.cam_select = s.cam_idx( boolFcn(FSArray) );
                s.epics_select = s.epics_idx( boolFcn(FSArray) );
                
                [~,idx] = intersect(s.cam_sort, s.cam_select);
                s.cam_idx = s.cam_sort(sort(idx));
                
                [~,idx2] = intersect(s.epics_sort, s.epics_select);
                s.epics_idx = s.epics_sort(sort(idx2));
                
            elseif nargin == 1
                s.cam_select = s.CAM_INDEX;
                s.epics_select = s.EPICS_INDEX;
                
                [~,idx] = intersect(s.cam_sort, s.cam_select);
                s.cam_idx = s.cam_sort(sort(idx));
                
                [~,idx2] = intersect(s.epics_select, s.epics_sort);
                s.epics_idx = s.epics_sort(sort(idx2));
            end
            
        end
        
        function visImages(s, diag)
            
            [data,diagData] = s.hlpCheckImage(s.dataSet, s.hdr, diag);

            %Pre-allocate
            figure(1)
            imagesc(diagData)
            w = waitforbuttonpress;
            

            for k = 2:s.visImageIncrement:length(data.dat)
                diagData = imread( sprintf('%s%s',s.hdr,data.dat{k}) );
                imagesc(diagData)
                w = waitforbuttonpress;
            end
            
        end
        
        function waterfallPlot(s, diag, fcn)

            if nargin < 3
                print('Need at least 3 inputs')
                return
            end

            [data,diagData] = s.hlpCheckImage(s.dataSet, s.hdr, diag);
            
            %Check that the provided function maps 2D -> 1D array
            wFData = fcn(diagData);

            [r,c] = size(wFData);
            if ~xor(r == 1, c == 1)
                error('The provided function does not map the data correctly. Mapped data has %d rows and %d columns.',r,d)
            end

            %Pre-allocate
            len = length(data.dat(s.cam_idx));
            waterfall = zeros(r*c,len);
            waterfall(:,1) = wFData;

            for k = 2:len
                diagData = imread( sprintf('%s%s',s.hdr,data.dat{s.cam_idx(k)}) );
                wFData = fcn(diagData);
                waterfall(:,k) = wFData;
            end

            figure;
            imagesc(waterfall)
            xlabel(s.sortLab)
            ylabel(func2str(fcn))
            %set(gca,'interpreter','none','fontsize',18)

        end
        
        function correlationPlot(s, inputArg1, inputArg2)
        %CORRELATIONPLOT plots 1 FACET-vector as a function another
        %  A FACET-vector can be defined by:
        %   - experimental scalar value captured for each shot in a dataset
        %   - a combination of an image diagnostic and a function that maps images
        %   to scalar values      
        % 
        %  Example use:
        %  
        %
        %  A scalar can be either a cell array containing a string with the a
        %  string with the name of a scalar diagnostic, or a cell containing the
        %  name of an image diagnostic and a function that maps the image to a
        %  scalar value.
        %  
        %  Single scalar case:
        %
        %  correlationPlot({'BPMS_LI20_2445_X'})
        %
        %    plots the scalar diagnostic value BPMS_LI20_2445_X as a 
        %  function of shot number
        %
        %  correlationPlot({'EOS_LO',@(x) sum(sum(x)) })
        %    plots the scalar diagnostic value from the image diagnostic 'EOS_LO' 
        %  mapped to a scalar by the lambda function @(x) sum(sum(x)) as a function 
        %  of shot number
        %  
        %  Double scalar case:
        %    correlationPlot({'BPMS_LI20_2445_X'}, {'EOS_LO',@(x) sum(sum(x))}) 
        %  plots the scalar diagnostic value BPMS_LI20_2445_X as a
        %  function of the scalar diagnostic value from the image diagnostic
        %  'EOS_LO' mapped to a scalar by the lambda function @(x) sum(sum(x))
        %
        %  Author: Henrik Ekerfelt
        %

        %% Input parsing
        
            if ischar(inputArg1)
                inputArg1 = {inputArg1};
            end
            type = s.hlpIsFSArray(inputArg1);
            
            
            if type
                [x, xlab] = s.hlpGetScalarArray(inputArg1,type);
            else 
                error('input argument #1 not a FACET Scalar Array')
            end
            
            if nargin == 2
                figure
                plot(x,'x')
                ylabel(xlab)
                xlabel('idx')
                %set(gca,'interpreter','none','fontsize',18)
                return
            end
            
            if nargin == 3
                type = s.hlpIsFSArray(inputArg2);
                if type
                    [y, ylab] = s.hlpGetScalarArray(inputArg2,type);
                else 
                    error('input argument #2 not a FACET Scalar Array')
                end
                figure
                plot(x,y,'x')
                ylabel(xlab)
                xlabel(ylab)
                %set(gca,'interpreter','none','fontsize',18)
                return
            end

end
        
        
    end
    
    methods (Access = private)
        
        function select_cams(s)
        %% Checks the UID of the images and the scalar data
        % selects and stores the matching indices in s.CAM_INDEX and
        % s.EPICS_INDEX. This is performeds once per initialization.
            
            CAMS = fields(s.dataSet.images);
            nCAMS = numel(CAMS);

            EPICS_UID  = s.dataSet.scalars.PATT_SYS1_1_PULSEID.UID;
            COMMON_UID = EPICS_UID;

            for i = 1:nCAMS
                cam_struct = s.dataSet.images.(CAMS{i});
                COMMON_UID = intersect(COMMON_UID,cam_struct.UID);
            end

            [~,~,EPICS_INDEX] = intersect(COMMON_UID,EPICS_UID);

            MATCHED     = numel(EPICS_INDEX);
            CAM_INDEX   = zeros(MATCHED,nCAMS);

            for i = 1:nCAMS
                cam_struct = s.dataSet.images.(CAMS{i});
                [~,~,CAM_INDEX(:,i)] = intersect(COMMON_UID,cam_struct.UID);
            end

            n_req   = s.dataSet.metadata.param.n_step*s.dataSet.metadata.param.n_shot;
            n_UID   = MATCHED;
            percent = 100*n_UID/n_req;
            per_str = num2str(percent,'%2.1f');
            
            s.EPICS_INDEX = EPICS_INDEX;
            if isequal(CAM_INDEX(:,1), unique(CAM_INDEX))
                s.CAM_INDEX = unique(CAM_INDEX);
            else
                warning('not all image diagnostics at all shots were saved, not handled well')
            end
            
            s.cam_idx = s.CAM_INDEX;
            s.epics_idx = s.EPICS_INDEX;
            s.cam_select = s.CAM_INDEX;
            s.epics_select = s.EPICS_INDEX;
            s.cam_sort = s.CAM_INDEX;
            s.epics_sort = s.EPICS_INDEX;
            s.include_data = 1:length(s.epics_idx)

            display([per_str '% of shots remain after UID matching']);
            
        end
        
        function type = hlpIsFSArray(s, inputArg)
        %% Determines if inputArg is a FACET Scalar Array and if so, what 
        % type.
        %
        % 0 - > not a FACET Scalar Array
        % 1 - > name of scalar diagnostic
        % 2 - > name of image diagnostic and a 2D - > scalar function
        %

            if ischar(inputArg)
                inputArg = {inputArg};
            end
            
            if ~iscell(inputArg)
                type = 0;
                %disp('not a cell in hlpIsFSArray')
                return
            end
            
            % Type 1 check
            if length(inputArg) == 1 && ischar(inputArg{1}) && ...
                    isfield(s.dataSet.scalars, inputArg{1})
                type = 1;
                return
            end
            
            % Type 2 check
            if length(inputArg) == 2 && ischar(inputArg{1}) && ...
                    isfield(s.dataSet.images, inputArg{1}) && ...
                    isa(inputArg{2},'function_handle') && ...
                    isscalar(inputArg{2}([1,2;3,4]))
                
                type = 2;
                return
            end
            
            type = 0;
            
        end
        
        function [scalarArray, scalarLabel] = hlpGetScalarArray(s, FACETscalar, type)
        %% hlpGetScalarArray extracts a 1D array of values from FACET data
        %  currently supports scalar diagnostics (type 1) and image
        %  diagnostics combined with a function that maps 2D - > scalar
        %  (type 2).

            if type == 1
                scalarArray = getfield(s.dataSet.scalars,FACETscalar{1}).dat(s.epics_idx);
                scalarLabel = FACETscalar{1};
                return
            elseif type == 2
                dS = getfield(s.dataSet.images,FACETscalar{1});
                nbrOfShots = length(s.cam_idx);
                scalarArray = zeros(1,nbrOfShots);
                fcn = FACETscalar{2};

                for k = 1:nbrOfShots
                    diagData = imread( sprintf('%s%s',s.hdr,dS.dat{s.cam_idx(k)}) );
                    scalarData = fcn(diagData);
                    scalarArray(:,k) = scalarData;
                end
                
                scalarLabel = sprintf('%s|%s',FACETscalar{1},func2str(fcn));
            end

        end
        
        function [data, diagData] = hlpCheckImage(s, dS,hdr,diag)


            % Check that diagnostic exists
            if ~isfield(dS.images,diag)
                error(sprintf('Could not find %s as an image diagnostic.',diag))
            end


            % Check that image links exist
            data = getfield(dS.images,diag);
            if isempty(data.dat)
                error(sprintf('In %s, images.%s.dat is empty',diag,diag))
            end

            %Find file format
            expr = '\.[0-9a-z]+$';
            idx = regexp(data.dat{1},expr) + 1;

            %Check and select method to load from file format
            if strcmp(data.dat{1}(idx:end),'tif') || strcmp(data.dat{1}(idx:end),'tiff')
                diagData = imread( sprintf('%s%s',hdr,data.dat{1}) );
            else
                error('File format not recognized for data.')
            end

        end

    end
end

