
%% ================================================================
%  PS-DSCI: Adaptive Detection and Noise-Corrected BFI Reconstruction
%
%  Features:
%    1. Adaptive illumination-source detection
%    2. Signed Euclidean distance mapping
%    3. SDS-dependent adaptive detector selection
%    4. Spatial speckle contrast calculation
%    5. Intensity-weighted noise correction
%    6. BFI reconstruction at multiple SDS values
%    7. Publication-quality visualization
%
%  Required:
%    MATLAB R2022b or later
%    Image Processing Toolbox
%
%  Input:
%    Image_(1).tiff, Image_(2).tiff, ...
%    dark_image.mat containing mean_img and std_img
%
%  Output:
%    Noise-corrected BFI maps (.mat)
%    Optional diagnostic figures
%    Optional source/detector boundary figures
%    Multi-SDS flow-map figures
%
%  NOTE:
%    BFI = 1/K_corr^2 is a relative flow index, not absolute
%    cerebral blood flow.
%% ================================================================

clc;
close all;
clear;

%% ================================================================
%  1. USER-DEFINED PATHS
%% ================================================================

% Set the folder containing your TIFF images
rootDirectory = 'Input Directory';

% Dark-image file is located one folder above rootDirectory
darkPath = fullfile(fileparts(rootDirectory), ...
                    'dark_image.mat');

% Results directory
resultsDir = fullfile(rootDirectory, 'Results');

if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

%% ================================================================
%  2. RECONSTRUCTION PARAMETERS
%% ================================================================

% Image numbers to process
frameList = 1:100;

% Source-detector separations (mm)
SD_list = [1,3,5];

% Detector-band thickness (mm)
thickness_list = 4;

% Speckle contrast window
windowSize = 5;

% Spatial calibration (pixels/mm)
radiusConversionF = 17;

% Detector-band overlap (pixels)
overlap = 20;

% Preserve original frame-skipping behavior
skipEvery100thFrame = true;

%% ================================================================
%  3. NOISE CORRECTION PARAMETERS
%% ================================================================

applyNoiseCorrection = true;

gain_db = 0;
gain_linear = 10^(gain_db/20);
gain = gain_linear * 0.0236;

%% ================================================================
%  4. SOURCE DETECTION PARAMETERS
%% ================================================================

gaussianSigma = 2;

thresholdFactor = 1.2;

erosionRadius = 5;

morphRadius = 5;

minimumObjectArea = 10;

minimumSourceArea = 100;

sourceSmoothRadius = 5;

%% ================================================================
%  5. FIGURE DISPLAY AND SAVE OPTIONS
%% ================================================================

% ------------------------------------------------
% Diagnostic figure:
% Raw intensity / source / distance / detector / K
% ------------------------------------------------

showDiagnosticFigures = false;
saveDiagnosticFigures = false;

% Save the five individual panels separately
saveIndividualPanels = false;

% ------------------------------------------------
% Source and detector boundary visualization
% ------------------------------------------------

showBoundaryFigures = false;
saveBoundaryFigures = false;

% SDS values for detector-boundary visualization
boundarySDList = [1 3 5];

% ------------------------------------------------
% Multi-SDS flow maps
% ------------------------------------------------

showFlowFigures = true;
saveFlowFigures = true;

% Save reconstructed data
saveFlowData = true;

% ------------------------------------------------
% Frames used for diagnostic visualization
%
% Flow reconstruction uses ALL frames in frameList.
% These settings only control intermediate figures.
% ------------------------------------------------

diagnosticFrames = 1;

boundaryFrames = 1;

%% ================================================================
%  6. FLOW VISUALIZATION SETTINGS
%% ================================================================

% Normalize all maps against ONE common reference
normalizeFlow = false;

% Normalized flow color limits
flowCLimNormalized = [0 2];

% Raw BFI limits when normalization is disabled
% [] means automatic shared limits
flowCLimRaw = [5 95];

% Publication figure settings
fontName = 'Arial';

titleFontSize = 18;
axisFontSize = 15;
colorbarFontSize = 16;

figureWidth = 1800;
figureHeight = 550;

exportResolution = 600;

% Optional scale bar
showScaleBar = true;
scaleBarMM = 5;

% Image cropping
useCrop = true;

cropRows = 300:1100;
cropCols = 350:1500;

% Diagnostic color limits
rawCLim = [0 256];
distanceCLim = [-10 10];
contrastCLim = [0 0.3];

%% ================================================================
%  7. LOAD DARK IMAGE
%% ================================================================

if ~exist(darkPath, 'file')
    error('Dark image not found: %s', darkPath);
end

darkData = load(darkPath);

assert(isfield(darkData, 'mean_img') && ...
       isfield(darkData, 'std_img'), ...
       'Dark file must contain mean_img and std_img.');

% Preserve the original 960 x 1256 dark-image region
mean_dark = double(darkData.mean_img);
std_dark  = double(darkData.std_img);

%% ================================================================
%  8. INITIALIZATION
%% ================================================================

nSD = numel(SD_list);
nThickness = numel(thickness_list);

% Accumulation arrays for each SDS and thickness
sumK = cell(nThickness, nSD);
X_mask = cell(nThickness, nSD);

% Final reconstructed maps
flowMaps = cell(nThickness, nSD);
contrastMaps = cell(nThickness, nSD);

fprintf('\n============================================\n');
fprintf('PS-DSCI Reconstruction\n');
fprintf('============================================\n');

fprintf('Input: %s\n', rootDirectory);

%% ================================================================
%  9. FRAME PROCESSING
%
%  Source detection and speckle contrast are calculated
%  once per frame and reused for all SDS values.
%% ================================================================

for counter = frameList

    if skipEvery100thFrame && mod(counter,100) == 0
        continue;
    end

    fileName = fullfile(rootDirectory, ...
        sprintf('Image_(%d).tiff', counter));

    if ~exist(fileName, 'file')
        warning('Missing frame %d. Skipping.', counter);
        continue;
    end

    %% ------------------------------------------------------------
    %  9.1 Load image
    %% ------------------------------------------------------------

    image_line = imread(fileName);

    image_raw = double(image_line);

    if ~isequal(size(image_raw), size(mean_dark))
        error(['Image and dark-image dimensions differ. ' ...
               'Check dark-image cropping.']);
    end

    % Subtract dark offset
    I = image_raw - mean_dark;

    [nRows, nCols] = size(I);

    %% ------------------------------------------------------------
    %  9.2 Initialize accumulation arrays
    %% ------------------------------------------------------------

    if isempty(sumK{1,1})

        for t = 1:nThickness
            for s = 1:nSD

                sumK{t,s} = zeros(nRows,nCols);

                X_mask{t,s} = zeros(nRows,nCols);

            end
        end
    end

    %% ------------------------------------------------------------
    %  9.3 Adaptive source extraction
    %% ------------------------------------------------------------

    img = imgaussfilt(image_line, gaussianSigma);

    level = graythresh(img);

    level = min(level * thresholdFactor, 1);

    image_binary = imbinarize(img, level);

    image_binary = imclose(image_binary, ...
        strel('disk', morphRadius));

    image_binary = imfill(image_binary, 'holes');

    image_binary = bwareaopen(image_binary, ...
        minimumObjectArea);

    labeledImage = bwlabel(image_binary, 8);

    stats = regionprops(labeledImage, 'Area');

    if isempty(stats)
        warning('No source detected in frame %d.', counter);
        continue;
    end

    allAreas = [stats.Area];

    validIdx = find(allAreas > minimumSourceArea);

    if isempty(validIdx)
        warning('No valid source in frame %d.', counter);
        continue;
    end

    source = ismember(labeledImage, validIdx);

    % Smooth source boundary
    source = imopen(source, ...
        strel('disk', sourceSmoothRadius));

    source = imclose(source, ...
        strel('disk', sourceSmoothRadius));

    source = imfill(source, 'holes');

    % Original source erosion
    source = imerode(source, ...
        strel('disk', erosionRadius));

    if ~any(source(:))
        warning('Empty source after erosion, frame %d.',counter);
        continue;
    end

    %% ------------------------------------------------------------
    %  9.4 Signed Euclidean distance map
    %% ------------------------------------------------------------

    D_out = bwdist(source);

    D_in = bwdist(~source);

    D_signed = D_out;

    D_signed(source) = -D_in(source);

    D_signed_mm = D_signed / radiusConversionF;

    %% ------------------------------------------------------------
    %  9.5 Noise-corrected speckle contrast
    %
    %  Calculated ONCE per frame
    %% ------------------------------------------------------------

    [K, meanI] = movingContrast( ...
        I, windowSize, std_dark, gain, ...
        applyNoiseCorrection);

    %% ------------------------------------------------------------
    %  9.6 Determine visualization crop
    %% ------------------------------------------------------------

    if useCrop

        % Clamp requested crop to available image dimensions
        rows = cropRows(cropRows >= 1 & cropRows <= nRows);

        cols = cropCols(cropCols >= 1 & cropCols <= nCols);

        if isempty(rows) || isempty(cols)
            error('Crop region is outside the image.');
        end

    else

        rows = 1:nRows;
        cols = 1:nCols;

    end

    %% ------------------------------------------------------------
    %  9.7 Reconstruct detector masks for each SDS
    %% ------------------------------------------------------------

    for s = 1:nSD

        SD_separation = SD_list(s);

        SD_pixels = SD_separation * radiusConversionF;

        for t = 1:nThickness

            thickness = thickness_list(t);

            belt = thickness * radiusConversionF;

            % Preserve original detector-band definition
            lowerBound = SD_pixels - overlap;

            upperBound = SD_pixels - overlap + belt;

            mask = (D_signed >= lowerBound) & ...
                   (D_signed <= upperBound);

            % Accumulate corrected contrast
            sumK{t,s} = sumK{t,s} + K .* double(mask);

            X_mask{t,s} = X_mask{t,s} + double(mask);

            %% ----------------------------------------------------
            %  Optional five-panel diagnostic figure
            %% ----------------------------------------------------

            if (showDiagnosticFigures || ...
                    saveDiagnosticFigures || ...
                    saveIndividualPanels) && ...
                    ismember(counter, diagnosticFrames)

                fig = figure( ...
                    'Color','w', ...
                    'Visible',visibilityString(showDiagnosticFigures), ...
                    'Position',[100 100 2100 480]);

                tl = tiledlayout(fig,1,5, ...
                    'TileSpacing','compact', ...
                    'Padding','compact');

                rawDisp = image_raw(rows,cols);
                sourceDisp = source(rows,cols);
                distDisp = D_signed_mm(rows,cols);
                maskDisp = mask(rows,cols);

                KDisp = K(rows,cols);
                KDisp(~maskDisp) = 0;

                panels = {rawDisp, sourceDisp, ...
                          distDisp, maskDisp, KDisp};

                panelTitles = { ...
                    '(a) Raw intensity', ...
                    '(b) Detected source', ...
                    '(c) Signed distance (mm)', ...
                    '(d) Detector mask', ...
                    '(e) Corrected contrast'};

                panelCLim = { ...
                    rawCLim, [0 1], distanceCLim, ...
                    [0 1], contrastCLim};

                panelMaps = { ...
                    gray(256), gray(256), parula(256), ...
                    gray(256), hot(256)};

                axList = gobjects(1,5);

                for p = 1:5

                    ax = nexttile(tl,p);
                    axList(p) = ax;

                    imagesc(ax,panels{p});

                    axis(ax,'image','off');

                    colormap(ax,panelMaps{p});

                    clim(ax,panelCLim{p});

                    title(ax,panelTitles{p}, ...
                        'FontName',fontName, ...
                        'FontSize',titleFontSize, ...
                        'FontWeight','bold');

                    cb = colorbar(ax);

                    styleColorbar(cb,colorbarFontSize);

                    if p == 2 || p == 4
                        cb.Ticks = [0 1];
                    end

                end

                prefix = sprintf( ...
                    'Frame_%d_SDS_%g_Thickness_%g', ...
                    counter,SD_separation,thickness);

                if saveDiagnosticFigures

                    exportgraphics(fig, ...
                        fullfile(resultsDir, ...
                        [prefix '_Diagnostic.png']), ...
                        'Resolution',exportResolution);

                end

                if saveIndividualPanels

                    names = {'Raw','Source','Distance', ...
                             'Detector','Contrast'};

                    for p = 1:5

                        exportgraphics(axList(p), ...
                            fullfile(resultsDir, ...
                            sprintf('%s_%s.png', ...
                            prefix,names{p})), ...
                            'Resolution',exportResolution);

                    end

                end

                if ~showDiagnosticFigures
                    close(fig);
                end

            end

        end
    end

    %% ------------------------------------------------------------
    %  9.8 Source and detector boundary visualization
    %% ------------------------------------------------------------

    if (showBoundaryFigures || saveBoundaryFigures) && ...
            ismember(counter,boundaryFrames)

        for t = 1:nThickness

            thickness = thickness_list(t);

            belt = thickness * radiusConversionF;

            nBoundary = numel(boundarySDList);

            figB = figure( ...
                'Color','w', ...
                'Visible',visibilityString(showBoundaryFigures), ...
                'Position',[100 100 400*nBoundary 520]);

            tlB = tiledlayout(figB,1,nBoundary, ...
                'TileSpacing','compact', ...
                'Padding','compact');

            for j = 1:nBoundary

                currentSDS = boundarySDList(j);

                currentSD = currentSDS * radiusConversionF;

                lowerBound = currentSD - overlap;

                upperBound = currentSD - overlap + belt;

                ax = nexttile(tlB,j);

                imagesc(ax,image_raw(rows,cols));

                colormap(ax,gray(256));

                clim(ax,rawCLim);

                axis(ax,'image','off');

                hold(ax,'on');

                % Source boundary (red)
                contour(ax,double(source(rows,cols)), ...
                    [0.5 0.5], ...
                    'Color',[1 0 0], ...
                    'LineWidth',2);

                % Inner detector boundary (green)
                contour(ax,D_signed(rows,cols), ...
                    [lowerBound lowerBound], ...
                    'Color',[0 0.8 0.2], ...
                    'LineWidth',2);

                % Outer detector boundary (cyan)
                contour(ax,D_signed(rows,cols), ...
                    [upperBound upperBound], ...
                    'Color',[0 0.8 1], ...
                    'LineWidth',2);

                title(ax, ...
                    sprintf('(%s) SDS = %g mm', ...
                    char('a'+j-1),currentSDS), ...
                    'FontName',fontName, ...
                    'FontSize',titleFontSize, ...
                    'FontWeight','bold');

                hold(ax,'off');

            end

            if saveBoundaryFigures

                exportgraphics(figB, ...
                    fullfile(resultsDir, ...
                    sprintf('Boundaries_Frame_%d_T_%g.png', ...
                    counter,thickness)), ...
                    'Resolution',exportResolution);

            end

            if ~showBoundaryFigures
                close(figB);
            end

        end
    end

    fprintf('Processed frame %d\n',counter);

end

%% ================================================================
%  10. FINAL BFI RECONSTRUCTION
%% ================================================================

for t = 1:nThickness

    for s = 1:nSD

        if isempty(sumK{t,s})
            warning('No valid frames for SDS=%g, T=%g.', ...
                SD_list(s),thickness_list(t));
            continue;
        end

        % Number of contributing detector samples
        counts = X_mask{t,s};

        validPixels = counts > 0;

        % Average corrected speckle contrast
        Kmean_2D = nan(size(counts));

        Kmean_2D(validPixels) = ...
            sumK{t,s}(validPixels) ./ counts(validPixels);

        % Reconstruct BFI
        flow_main = nan(size(Kmean_2D));

        positiveK = validPixels & ...
                    isfinite(Kmean_2D) & ...
                    (Kmean_2D > 0);

        flow_main(positiveK) = ...
            1 ./ (Kmean_2D(positiveK).^2);

        flowMaps{t,s} = flow_main;

        contrastMaps{t,s} = Kmean_2D;

        %% --------------------------------------------------------
        %  Save reconstructed data
        %% --------------------------------------------------------

        if saveFlowData

            SD_separation = SD_list(s);
            thickness = thickness_list(t);

            X_mask_current = counts;

            save(fullfile(resultsDir, ...
                sprintf('NC_SDS_%g_Thickness_%g.mat', ...
                SD_separation,thickness)), ...
                'flow_main','Kmean_2D', ...
                'X_mask_current', ...
                'SD_separation','thickness');

        end

    end
end

%% ================================================================
%  11. MULTI-SDS FLOW MAP VISUALIZATION
%
%  One figure per detector-band thickness
%  Horizontal subplots for different SDS values
%  Shared normalization reference and hot colorbar
%% ================================================================

if showFlowFigures || saveFlowFigures

    for t = 1:nThickness

        %% --------------------------------------------------------
        %  Determine common normalization reference
        %% --------------------------------------------------------

        % Use first valid SDS map as the reference
        referenceValue = NaN;

        for s = 1:nSD

            flow = flowMaps{t,s};

            if isempty(flow)
                continue;
            end

            flowCrop = flow(rows,cols);

            validValues = flowCrop( ...
                isfinite(flowCrop) & flowCrop > 0);

            if ~isempty(validValues)

                referenceValue = mean(validValues(:));
                break;

            end

        end

        if isnan(referenceValue)
            warning('No valid flow maps for thickness %g.', ...
                thickness_list(t));
            continue;
        end

        %% --------------------------------------------------------
        %  Determine common color limits
        %% --------------------------------------------------------

        if normalizeFlow

            displayCLim = flowCLimNormalized;
            colorbarLabel = 'Normalized BFI';

        else

            colorbarLabel = 'BFI (a.u.)';

            if isempty(flowCLimRaw)

                allValid = [];

                for s = 1:nSD

                    flow = flowMaps{t,s};

                    if isempty(flow)
                        continue;
                    end

                    flowCrop = flow(rows,cols);

                    values = flowCrop( ...
                        isfinite(flowCrop) & flowCrop > 0);

                    allValid = [allValid; values(:)];

                end

                if isempty(allValid)
                    continue;
                end

                displayCLim = prctile(allValid,[2 98]);

                if displayCLim(1) == displayCLim(2)
                    displayCLim = [0 max(allValid)];
                end

            else

                displayCLim = flowCLimRaw;

            end

        end

        %% --------------------------------------------------------
        %  Create publication-quality flow figure
        %% --------------------------------------------------------

        figFlow = figure( ...
            'Color','w', ...
            'Visible',visibilityString(showFlowFigures), ...
            'Position',[100 100 figureWidth figureHeight]);

        tlFlow = tiledlayout(figFlow,1,nSD, ...
            'TileSpacing','compact', ...
            'Padding','compact');

        axFlow = gobjects(1,nSD);

        for s = 1:nSD

            ax = nexttile(tlFlow,s);

            axFlow(s) = ax;

            flow = flowMaps{t,s};

            if isempty(flow)
                axis(ax,'off');
                continue;
            end

            flowCrop = flow(rows,cols);

            % Apply ONE common reference to all SDS maps
            if normalizeFlow
                flowDisplay = flowCrop / referenceValue;
            else
                flowDisplay = flowCrop;
            end

            imagesc(ax,flowDisplay);

            colormap(ax,hot(256));

            clim(ax,displayCLim);

            axis(ax,'image','off');

            title(ax, ...
                sprintf('(%s) SDS = %g mm', ...
                char('a'+s-1),SD_list(s)), ...
                'FontName',fontName, ...
                'FontSize',titleFontSize, ...
                'FontWeight','bold');

            %% ----------------------------------------------------
            %  Scale bar
            %% ----------------------------------------------------

            if showScaleBar

                scalePixels = scaleBarMM * radiusConversionF;

                xStart = 0.06 * numel(cols);

                yPos = 0.92 * numel(rows);

                hold(ax,'on');

                plot(ax, ...
                    [xStart xStart+scalePixels], ...
                    [yPos yPos], ...
                    'w-', ...
                    'LineWidth',4);

                text(ax, ...
                    xStart+scalePixels/2, ...
                    yPos-12, ...
                    sprintf('%g mm',scaleBarMM), ...
                    'Color','w', ...
                    'FontName',fontName, ...
                    'FontSize',14, ...
                    'FontWeight','bold', ...
                    'HorizontalAlignment','center');

                hold(ax,'off');

            end

        end

        %% --------------------------------------------------------
        %  Shared colorbar
        %% --------------------------------------------------------

        cb = colorbar(axFlow(end));

        cb.Layout.Tile = 'east';

        styleColorbar(cb,colorbarFontSize);

        cb.Label.String = colorbarLabel;

        cb.Label.FontSize = colorbarFontSize;

        cb.Label.FontWeight = 'bold';

        %% --------------------------------------------------------
        %  Save figure
        %% --------------------------------------------------------

        if saveFlowFigures

            thickness = thickness_list(t);

            exportgraphics(figFlow, ...
                fullfile(resultsDir, ...
                sprintf('Multi_SDS_Flow_T_%g.png',thickness)), ...
                'Resolution',exportResolution);

            % Optional vector PDF
            exportgraphics(figFlow, ...
                fullfile(resultsDir, ...
                sprintf('Multi_SDS_Flow_T_%g.pdf',thickness)), ...
                'ContentType','vector');

        end

        if ~showFlowFigures
            close(figFlow);
        end

    end
end

fprintf('\n============================================\n');
fprintf('Reconstruction complete.\n');
fprintf('Results: %s\n',resultsDir);
fprintf('============================================\n');

%% ================================================================
%  LOCAL FUNCTION 1:
%  Spatial speckle contrast and noise correction
%% ================================================================

function [K_corr,meanI] = movingContrast( ...
    I,windowSize,std_dark,gain,applyNoiseCorrection)

    kernel = ones(windowSize)/windowSize^2;

    meanI = conv2(I,kernel,'same');

    meanI2 = conv2(I.^2,kernel,'same');

    varI = max(meanI2-meanI.^2,0);

    % Raw speckle contrast
    K = sqrt(varI)./meanI;

    K(meanI <= 0) = 0;

    K(~isfinite(K)) = 0;

    if ~applyNoiseCorrection

        K_corr = K;
        return;

    end

    %% ------------------------------------------------------------
    %  Intensity-weighted noise correction
    %
    %  Preserves the original user-defined formulation
    %% ------------------------------------------------------------

    globalMean = mean(meanI(:));

    % Read noise
    v_read = (std_dark.^2) ./ (meanI.^2);

    % Shot noise
    v_shot = (gain .* meanI) ./ (meanI.^2);

    % Quantization noise
    v_quant = (1/12) ./ (meanI.^2);

    % Weighted noise variance
    v_noise= (v_read + v_shot + v_quant);

    % Corrected speckle contrast
    K_corr = sqrt(max(K.^2-v_noise,0));

    % Invalid intensity values
   K_corr(meanI <= 0) = NaN;

   K_corr(~isfinite(K_corr)) = NaN;

end

%% ================================================================
%  LOCAL FUNCTION 2:
%  Figure visibility
%% ================================================================

function state = visibilityString(showFigure)

    if showFigure
        state = 'on';
    else
        state = 'off';
    end

end

%% ================================================================
%  LOCAL FUNCTION 3:
%  Publication-quality colorbar formatting
%% ================================================================

function styleColorbar(cb,fontSize)

    cb.FontName = 'Arial';

    cb.FontSize = fontSize;

    cb.FontWeight = 'bold';

    cb.LineWidth = 1.5;

    cb.TickDirection = 'out';

    cb.Color = 'k';

end
