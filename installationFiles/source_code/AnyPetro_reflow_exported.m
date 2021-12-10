classdef AnyPetro_reflow_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        AnyPetroUIFigure                matlab.ui.Figure
        GridLayout                      matlab.ui.container.GridLayout
        LeftPanel                       matlab.ui.container.Panel
        GridLayout3                     matlab.ui.container.GridLayout
        SelectdatacolumnDropDown        matlab.ui.control.DropDown
        SelectdatacolumnDropDownLabel   matlab.ui.control.Label
        SelectcolumnfordatauseDropDownLabel  matlab.ui.control.Label
        SelectcolumnfordatauseDropDown  matlab.ui.control.DropDown
        SelectdatatypecolumnDropDown    matlab.ui.control.DropDown
        SelectdatatypecolumnDropDownLabel  matlab.ui.control.Label
        UITable                         matlab.ui.control.Table
        SelectcolumnswithsamplespecificsListBoxLabel  matlab.ui.control.Label
        SelectcolumnswithsamplespecificsListBox  matlab.ui.control.ListBox
        SelectdataweightscolumnDropDownLabel  matlab.ui.control.Label
        SelectdataweightscolumnDropDown  matlab.ui.control.DropDown
        uselog10ofdataCheckBox          matlab.ui.control.CheckBox
        GeneratedatavectorsButton       matlab.ui.control.Button
        VieworeditdatatableButton       matlab.ui.control.Button
        LoaddatafileButton              matlab.ui.control.Button
        CenterPanel                     matlab.ui.container.Panel
        GridLayout2                     matlab.ui.container.GridLayout
        ExportresultButton              matlab.ui.control.Button
        ParametertransformationButtonGroup  matlab.ui.container.ButtonGroup
        RangeButton                     matlab.ui.control.RadioButton
        LogButton                       matlab.ui.control.RadioButton
        NoneButton                      matlab.ui.control.RadioButton
        StartinversionButton            matlab.ui.control.Button
        InversionParametersPanel        matlab.ui.container.Panel
        GridLayout6                     matlab.ui.container.GridLayout
        ToleranceEditField              matlab.ui.control.NumericEditField
        ToleranceEditFieldLabel         matlab.ui.control.Label
        PerturbationEditField           matlab.ui.control.NumericEditField
        PerturbationEditFieldLabel      matlab.ui.control.Label
        MaximumnumberofiterationsEditField  matlab.ui.control.NumericEditField
        MaximumnumberofiterationsEditFieldLabel  matlab.ui.control.Label
        RegularizationPanel             matlab.ui.container.Panel
        GridLayout5                     matlab.ui.container.GridLayout
        lambda2EditField                matlab.ui.control.NumericEditField
        lambda2EditFieldLabel           matlab.ui.control.Label
        lambda1EditField                matlab.ui.control.NumericEditField
        lambda1EditFieldLabel           matlab.ui.control.Label
        lambda0EditField                matlab.ui.control.NumericEditField
        lambda0EditFieldLabel           matlab.ui.control.Label
        FOLabel                         matlab.ui.control.Label
        SelectForwardOperatorButton     matlab.ui.control.Button
        UITable2                        matlab.ui.control.Table
        RightPanel                      matlab.ui.container.Panel
        GridLayout4                     matlab.ui.container.GridLayout
        ClearcommandlineButton          matlab.ui.control.Button
        CommandlineoutputTextAreaLabel  matlab.ui.control.Label
        CommandlineoutputTextArea       matlab.ui.control.TextArea
        UIAxes2                         matlab.ui.control.UIAxes
        UIAxes                          matlab.ui.control.UIAxes
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
        twoPanelWidth = 768;
    end

    
    properties (Access = private)
        dataTable % Description
        DialogApp
        sampleSpecifics
        dataColumn
        dataWeightColumn
        dobs
        err
        options
        NoModPar
        modelTable
        ModelParameterNames
        ModelParameterWeights
        FOname
        mod_erg
        mod_std
        Plotname
        data_calc_final
        dataFilePath
        FOfilePath
    end
    
    methods (Access = private)
        
        function myUpdateAppLayout(app,event)
            currentFigureWidth = app.AnyPetroUIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 3x1 grid
                app.GridLayout.RowHeight = {'1x', '1x', '1x'};
                app.GridLayout.ColumnWidth = {'1x'};
                app.CenterPanel.Layout.Row = 2;
                app.CenterPanel.Layout.Column = 1;
                app.LeftPanel.Layout.Row = 1;
                app.LeftPanel.Layout.Column = 1;
                app.RightPanel.Layout.Row = 3;
                app.RightPanel.Layout.Column = 1;
            elseif (currentFigureWidth > app.onePanelWidth && currentFigureWidth <= app.twoPanelWidth)
                % Change to a 2x2 grid
                app.GridLayout.RowHeight = {'1x', '1x'};
                app.GridLayout.ColumnWidth = {'1x', '1x'};
                app.CenterPanel.Layout.Row = 1;
                app.CenterPanel.Layout.Column = 2;
                app.LeftPanel.Layout.Row = 1;
                app.LeftPanel.Layout.Column = 1;
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = [1,2];
            else
                % Change to a 1x3 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {'1x', '1x', '1x'};
                app.LeftPanel.Layout.Row = 1;
                app.LeftPanel.Layout.Column = 1;
                app.CenterPanel.Layout.Row = 1;
                app.CenterPanel.Layout.Column = 2;
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 3;
            end
        end
        
        function DoInversion(app)
            %% DOINVERSION %%
            app.options.transformation = app.ParametertransformationButtonGroup.SelectedObject.Text;
            modData = app.UITable2.Data;
            mod_start = modData.('startingValue');
            app.options.ParaRange = [modData.('lowerBound') modData.('upperBound')];
            app.options.mod_ref = modData.('referenceValue');
            app.options.modelParameterWeights = modData.('weight');
            app.options.maxi = app.MaximumnumberofiterationsEditField.Value;
            app.options.pert = app.PerturbationEditField.Value;
            app.options.tol = app.ToleranceEditField.Value;
            app.options.DataLog = app.uselog10ofdataCheckBox.Value;
            
            if isfield(app.options,'transformation')
                if strcmp(app.options.transformation,'Log') == 1
                    % work with logarithmic values
                    mod_start = log(mod_start);
                    mod_ref   = log(app.options.mod_ref);
                elseif strcmp(app.options.transformation,'Range') == 1
                    % work with constrained parameter ranges
                    mod_start = log(mod_start               - app.options.ParaRange(:,1)) ...
                        -log(app.options.ParaRange(:,2) - mod_start             );
                    mod_ref = log( app.options.mod_ref         - app.options.ParaRange(:,1)) ...
                        -log(app.options.ParaRange(:,2) -  app.options.mod_ref );
                elseif strcmp(app.options.transformation,'None') == 1
                    mod_ref   = app.options.mod_ref;
                else
                    error('Unknown tranformation.')
                end
            end
            
            % initial regularization parameter
            lambda              = [ app.lambda0EditField.Value; ...
                app.lambda1EditField.Value; ...
                app.lambda2EditField.Value];
            
            % synthetic data for starting model
            d = ForwardOperator(app,mod_start,app.options);
            
            % initial data residual
            data_residual       = app.dobs - d;
            
            %% data weight matrix Cd
            Cd = diag(app.err);
            
            %% regularization matrix Cm
            n       = numel(mod_start);
            Cm      = diag(app.options.modelParameterWeights);
            
            if lambda(2) ~= 0 || lambda(3) ~= 0
                ApplyC1C2 = modData.applyC1C2;
                if sum(ApplyC1C2) ~= 0
                    
                    if lambda(2) ~= 0
                        % include first derivative in regularization
                        C1      = 0.5.*eye(n-1,n) + [zeros(n-1,1) -0.5.*eye(n-1,n-1)];
                        C1(ApplyC1C2==0,:) = [];
                        C1(:,ApplyC1C2==0) = 0;
                        
                    else
                        C1 = zeros(n);
                    end
                    
                    if lambda(3) ~= 0
                        % include 2nd derivative in regularization
                        C2      = 0.25*(-eye(n-2,n) + [zeros(n-2,1) 2*eye(n-2,n-1)]  + [zeros(n-2,2) -eye(n-2,n-2)]);
                        C2(ApplyC1C2==0,:) = [];
                        C2(:,ApplyC1C2==0) = 0;
                    else
                        C2 = zeros(n);
                    end
                    
                    Cm       = [ Cm; lambda(2)*C1; lambda(3)*C2];
                end
            end
            
            %% initial objective function
            datanorm            = 0.5 * norm(Cd*    data_residual    )^2;
            modelnorm           = 0.5 * norm(Cm*(mod_start - mod_ref))^2;
            
            gamma = 0.00005;
            if modelnorm ~=0 && lambda(1)==0
                lambda(1) = gamma*(datanorm/modelnorm);
                app.options.lambdaStore = lambda(1);
            elseif modelnorm ==0 && lambda(1)==0
                lambda(1) = 1e-4;
            end
            
            obj_func_current    = datanorm + lambda(1)*modelnorm;
            obj_func            = obj_func_current;
            
            app.mod_erg             = mod_start;
            
            app.CommandlineoutputTextArea.Value = [app.CommandlineoutputTextArea.Value;
                {' '};
                {'------------------'};
                {'Starting inversion'};
                {'------------------'};
                {datestr(now,'yyyy-mm-dd_HH-MM-SS')};
                {'------------------'};
                {'options:'};
                {evalc('disp(app.options)')};
                {sprintf('Starting Gauss-Newton inversion')};
                {sprintf('Data residual norm:\t%1.2e',norm(data_residual,2))};
                {sprintf('Calculated lambda:\t%1.2e',lambda(1))};
                {sprintf('Objective function:\t%1.6e\n',obj_func_current)}];
            drawnow;
            app.CommandlineoutputTextArea.scroll('bottom');

            for i = 1:app.options.maxi
                app.CommandlineoutputTextArea.Value = [app.CommandlineoutputTextArea.Value;
                    {sprintf('Iteration: %d',i)}];
                drawnow;
                app.CommandlineoutputTextArea.scroll('bottom');
                app.options.mod_erg = app.mod_erg;
                %% Jacobian
                J                   = GetJacobian(app);
                
                %% compute model update
                
                delta_m             = SolveSystemOfNormalEquations(app,J,data_residual,lambda,Cd,Cm,mod_ref);
                
                %% line search
                step                = 1;
                mod_erg_temp        = app.mod_erg + step*delta_m;
                d_temp              = ForwardOperator(app,mod_erg_temp,app.options);
                obj_func_temp       = 0.5*norm(Cd*(app.dobs-d_temp))^2 + (lambda(1)/2)*norm(Cm*(mod_erg_temp - mod_ref))^2;
                diff_obj_func       = obj_func_current(end) - obj_func_temp;
                directional_deriv   = J*delta_m;
                
                if isnan(diff_obj_func)
                    diff_obj_func = 0.01*directional_deriv(1);
                end
                
                while diff_obj_func < 0.25*step*directional_deriv(1)  || ~isreal(d_temp)
                    step            = 0.5*step;
                    if step < 1e-6
                        app.CommandlineoutputTextArea.Value = [app.CommandlineoutputTextArea.Value;
                            {sprintf('Breakdown in line search. Exit.\n')}];
                        drawnow;
                        app.CommandlineoutputTextArea.scroll('bottom');
                        break
                    end
                    mod_erg_temp    = app.mod_erg + step*delta_m;
                    d_temp          = ForwardOperator(app,mod_erg_temp,app.options);
                    obj_func_temp   = 0.5*norm(Cd*(app.dobs-d_temp))^2 + (lambda(1)/2)*norm(Cm*(mod_erg_temp-mod_ref))^2;
                    diff_obj_func   = obj_func_current(end) - obj_func_temp;
                end
                
                dNormTemp = 0.5*norm(Cd*(app.dobs-d_temp))^2;
                mNormTemp = (lambda(1)/2)*norm(Cm*(mod_erg_temp-mod_ref))^2;
                
                app.CommandlineoutputTextArea.Value = [app.CommandlineoutputTextArea.Value;
                    {sprintf('Line search: Step size\t%1.6e',step)};
                    {sprintf('Data norm:\t\t\t%1.6e',dNormTemp)};
                    {sprintf('Model norm:\t\t\t%1.6e',mNormTemp)};
                    {sprintf('Objective function:\t\t%1.6e\n',obj_func_temp)}];
                app.CommandlineoutputTextArea.scroll('bottom');
                
                %% checking tolerance
                rel_change_obj_func = diff_obj_func / obj_func_current(end);
                
                %% updating
                d                   = d_temp;
                app.mod_erg         = mod_erg_temp;
                data_residual       = app.dobs - d;
                obj_func_current    = obj_func_temp;
                obj_func            = [obj_func; obj_func_current];
                
                %% update plot 1
                semilogy(app.UIAxes,[0:(length(obj_func)-1)],obj_func./obj_func(1),'s-')
                set(app.UIAxes,'PlotBoxAspectRatio',[1 1 1])
                
                %% update plot 2
                app.UIAxes2.XLimMode = 'auto';
                app.UIAxes2.YLimMode = 'auto';
                scatter(app.UIAxes2,app.dobs,d_temp,20,app.options.DataType,'o')
                colormap(app.UIAxes2,lines)
                set(app.UIAxes2,'PlotBoxAspectRatio',[1 1 1])
                hold(app.UIAxes2,'on')
                
                if ~app.uselog10ofdataCheckBox.Value
                    set(app.UIAxes2,'XScale','log','YScale','log')
                end
                
                app.UIAxes2.XLim = [min([app.UIAxes2.XLim(1);app.UIAxes2.YLim(1)]) ...
                    max([app.UIAxes2.XLim(2);app.UIAxes2.YLim(2)])];
                app.UIAxes2.YLim = app.UIAxes2.XLim;
                plot(app.UIAxes2,app.UIAxes2.XLim,app.UIAxes2.XLim,'k--')
                hold(app.UIAxes2,'off')
                drawnow;
                
                if rel_change_obj_func < app.options.tol
                    
                    app.CommandlineoutputTextArea.Value = [app.CommandlineoutputTextArea.Value;
                        {sprintf('Relative misfit change %1.6e is smaller than %1.6e. Stopping.\n', ...
                        rel_change_obj_func,app.options.tol)}];
                    drawnow;
                    app.CommandlineoutputTextArea.scroll('bottom');
                    break
                elseif isnan(dNormTemp) || isnan(mNormTemp) || isnan(obj_func_current)
                    app.CommandlineoutputTextArea.Value = [app.CommandlineoutputTextArea.Value;
                        {sprintf('Something went wrong (detecting NaN values). Stopping.\n')}];
                    break
                end
            end
            
            if app.options.DataLog
                app.data_calc_final = 10.^(d);
            else
                app.data_calc_final = d;
            end
            
            % estimate model confidence (after Malecki et al. 2020, GJI)
            Cd = (Cd')*Cd;
            redundancy = trace( eye(size(Cd,1)) - J*((J'*Cd*J)^(-1))*J'*Cd );
            s20 = ( (J*(step*delta_m)-data_residual)' ...
                *Cd *(J*(step*delta_m)-data_residual) )  /redundancy;
            CofactorMatrix = (J'*Cd*J+ lambda(1)*(Cm'*Cm))^(-1);

            sii = CofactorMatrix( sub2ind(size(CofactorMatrix) , ...
                      1:size(CofactorMatrix,1), ...
                      1:size(CofactorMatrix,2)));

            empiricalVariance = s20*sii';
            
            % go back to normal parameters
            if isfield(app.options,'transformation')
                if strcmp(app.options.transformation,'Log') == 1
                    % work with logarithmic values
                    empiricalVariance = (exp(app.mod_erg).^2).*empiricalVariance;
                    app.mod_std = sqrt(empiricalVariance);
        
                    app.mod_erg = exp(app.mod_erg);
                elseif strcmp(app.options.transformation,'Range') == 1
                    % work with constrained parameter ranges
                    empiricalVariance = ...
                        ((app.options.ParaRange(:,2).*exp(app.mod_erg)) ...
                        ./(exp(app.mod_erg) + 1) ...
                        - ( exp(app.mod_erg).*(app.options.ParaRange(:,1) ...
                        + app.options.ParaRange(:,2).*exp(app.mod_erg)) ) ...
                        ./ (exp(app.mod_erg) + 1).^2).^2 .*empiricalVariance;
%                     app.mod_std = sqrt(empiricalVariance);
                    app.mod_std = sqrt(empiricalVariance);
                    
                    app.mod_erg = (exp(app.mod_erg).*app.options.ParaRange(:,2) ...
                        + app.options.ParaRange(:,1))./(1+exp(app.mod_erg));
                    
                end
            end
            
            app.CommandlineoutputTextArea.Value = [app.CommandlineoutputTextArea.Value;
                {sprintf('Data residual norm after %d iterations:\t%1.2e',i,norm(app.dobs-d,2))};
                {sprintf('Value of objective function:\t\t\t%1.6e',obj_func_current(end))};
                {sprintf('Lambda:\t\t\t\t\t\t\t%1.2e\n',lambda(1))};
                {sprintf('Inversion result:')}];
            
            for pl = 1:numel(app.mod_erg)
                app.CommandlineoutputTextArea.Value = [app.CommandlineoutputTextArea.Value;
                    {sprintf('%s\t= %.2e +/- %.2e',app.ModelParameterNames{pl},app.mod_erg(pl),app.mod_std(pl))}];
            end
            app.CommandlineoutputTextArea.Value = [app.CommandlineoutputTextArea.Value;
                {' '}
                {datestr(now,'yyyy-mm-dd_HH-MM-SS')}];
            drawnow;
            app.CommandlineoutputTextArea.scroll('bottom');
            
            %%
            app.ExportresultButton.Enable = 'on';
            uialert(app.AnyPetroUIFigure,'Inversion finished.','Inversion',"Icon","success");
            figure(app.AnyPetroUIFigure)
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function delta_m = SolveSystemOfNormalEquations(app,J,data_residual,lambda,Cd,Cm,mod_ref)
            %% SOLVE %%
            JTJ     = J'*(Cd')*Cd*J;
            
            A       = JTJ + lambda(1)*(Cm'*Cm);
            B       = J'*(Cd')*Cd*data_residual - lambda(1)*(Cm'*Cm)*(app.mod_erg -mod_ref);
            delta_m = A\B;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function J = GetJacobian(app)
            m       = app.mod_erg;
            da_0    = ForwardOperator(app,m,app.options);
            m_neu   = m;
            
            for i=1:length(m)
                if m(i) == 0
                    m_mod = app.options.pert;
                    d_m= app.options.pert;
                else
                    m_mod = (1+ app.options.pert)*m(i);
                    d_m= app.options.pert*m(i);
                end
                m_neu(i) = m_mod;
                
                da_neu = ForwardOperator(app,m_neu,app.options);
                
                J(:,i) = (da_neu-da_0)/(d_m);
                
                m_neu(i) = m(i);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function dcalc = ForwardOperator(app,mod,options)
            if isfield(options,'transformation')
                if strcmp(options.transformation,'Log') == 1
                    mod = exp(mod);
                elseif strcmp(options.transformation,'Range') == 1
                    mod = (exp(mod).*app.modelTable.upperBound + app.modelTable.lowerBound)./(1+exp(mod));
                end
            end
            
            if isfield(options,'AuxiliaryStatements')
                for i = 1:numel(options.AuxiliaryStatements)
                    eval(options.AuxiliaryStatements{i});
                end
            end
            
            dcalc = nan(numel(options.DataType),1);
            
            
            for i = 1:numel(dcalc)
                if isfield(options,'FOinput') || isfield(options,'AddInput')
                    if isfield(options,'FOinput')
                        FOinput = options.FOinput(options,i)';
                        FOinput = num2cell(FOinput);
                        
                        if isfield(options,'AddInput')
                            AddInput = eval(options.AddInput{:});
                            ExtraInput = [mat2cell(AddInput,max(size(AddInput))), ... 
                                FOinput];
                        else
                            ExtraInput = FOinput;
                        end
                    else
                        AddInput = eval(options.AddInput{:});
                        ExtraInput = mat2cell(AddInput,ones(numel(AddInput),1));
                    end
                    dcalc(i) = options.FO{options.DataType(i)}(mod,ExtraInput{:});
                else
                    dcalc(i) = options.FO{options.DataType(i)}(mod);
                end
            end
            
            if options.DataLog
                dcalc = log10(dcalc);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    methods (Access = public)
        
        function updateDataTable(app,updatedDataTable)
            app.dataTable = updatedDataTable;
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            clc
            app.AnyPetroUIFigure.SizeChangedFcn = createCallbackFcn(app, @myUpdateAppLayout, true);
            app.GridLayout.ColumnWidth = {'1x', '1x', '1x'};
            app.lambda0EditFieldLabel.Text = [char(955) '_0'];
            app.lambda1EditFieldLabel.Text = [char(955) '_1'];
            app.lambda2EditFieldLabel.Text = [char(955) '_2'];
            
        end

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.AnyPetroUIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 3x1 grid
                app.GridLayout.RowHeight = {610, 610, 610};
                app.GridLayout.ColumnWidth = {'1x'};
                app.CenterPanel.Layout.Row = 1;
                app.CenterPanel.Layout.Column = 1;
                app.LeftPanel.Layout.Row = 2;
                app.LeftPanel.Layout.Column = 1;
                app.RightPanel.Layout.Row = 3;
                app.RightPanel.Layout.Column = 1;
            elseif (currentFigureWidth > app.onePanelWidth && currentFigureWidth <= app.twoPanelWidth)
                % Change to a 2x2 grid
                app.GridLayout.RowHeight = {610, 610};
                app.GridLayout.ColumnWidth = {'1x', '1x'};
                app.CenterPanel.Layout.Row = 1;
                app.CenterPanel.Layout.Column = [1,2];
                app.LeftPanel.Layout.Row = 2;
                app.LeftPanel.Layout.Column = 1;
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 2;
            else
                % Change to a 1x3 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {362, '1x', 443};
                app.LeftPanel.Layout.Row = 1;
                app.LeftPanel.Layout.Column = 1;
                app.CenterPanel.Layout.Row = 1;
                app.CenterPanel.Layout.Column = 2;
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 3;
            end
        end

        % Button pushed function: LoaddatafileButton
        function LoaddatafileButtonPushed(app, event)
            app.options = [];
            if isempty(app.dataFilePath)
                [dataFile,app.dataFilePath,~] = uigetfile({'*.xlsx';'*.txt';'*.mat';'*.*'},'Select data file');
            else
                [dataFile,app.dataFilePath,~] = uigetfile({'*.xlsx';'*.txt';'*.mat';'*.*'},'Select data file',app.dataFilePath);
            end
            
            if ~isequal(dataFile,0)
                [ ~, ~,ext] = fileparts([app.dataFilePath '\' dataFile]);
                if strcmp(ext,'.mat')
                    app.dataTable = load([app.dataFilePath '\' dataFile]);
                    app.dataTable = app.dataTable.(cell2mat(fieldnames(app.dataTable)));
                elseif strcmp(ext,'.xlsx')
                    app.dataTable = readtable([app.dataFilePath '\' dataFile]);
                elseif strcmp(ext,'.txt')
                    app.dataTable = readtable([app.dataFilePath '\' dataFile]);
                else
                    app.CommandlineoutputTextArea.Value = [app.CommandlineoutputTextArea.Value;
                        {sprintf('Unknown data file format. Import might fail.')}];
                    app.dataTable = readtable([app.dataFilePath '\' dataFile]);
                end
                
                app.VieworeditdatatableButton.Enable = true;
                
                numericColumns = varfun(@isnumeric,app.dataTable,'output','uniform');
                app.SelectcolumnswithsamplespecificsListBox.Items = sort(app.dataTable.Properties.VariableNames);
                app.SelectcolumnswithsamplespecificsListBoxLabel.Enable = true;
                app.SelectcolumnswithsamplespecificsListBox.Enable = true;
                app.SelectcolumnswithsamplespecificsListBox.Value = app.SelectcolumnswithsamplespecificsListBox.Items;
                app.sampleSpecifics = app.SelectcolumnswithsamplespecificsListBox.Value;
                
                app.SelectdatatypecolumnDropDown.Items = sort(app.dataTable.Properties.VariableNames);
                app.SelectdatatypecolumnDropDownLabel.Enable = true;
                app.SelectdatatypecolumnDropDown.Enable = true;
                typ = find(cellfun(@(x) ~isempty(x),strfind(lower(app.SelectdatatypecolumnDropDown.Items),'type')));
                if ~isempty(typ)
                    app.SelectdatatypecolumnDropDown.Value = app.SelectdatatypecolumnDropDown.Items{typ(1)};
                else
                    app.SelectdatatypecolumnDropDown.Value = app.SelectdatatypecolumnDropDown.Items{1};
                end
                
                app.SelectcolumnfordatauseDropDown.Items = sort(app.dataTable.Properties.VariableNames);
                app.SelectcolumnfordatauseDropDownLabel.Enable = true;
                app.SelectcolumnfordatauseDropDown.Enable = true;
                use = find(cellfun(@(x) ~isempty(x),strfind(lower(app.SelectcolumnfordatauseDropDown.Items),'use')));
                if ~isempty(use)
                    app.SelectcolumnfordatauseDropDown.Value = app.SelectcolumnfordatauseDropDown.Items{use(1)};
                else
                    app.SelectcolumnfordatauseDropDown.Value = app.SelectcolumnfordatauseDropDown.Items{1};
                end
                
                app.SelectdataweightscolumnDropDownLabel.Enable = true;
                app.SelectdataweightscolumnDropDown.Enable = true;
                app.SelectdataweightscolumnDropDown.Items = sort(app.dataTable.Properties.VariableNames(numericColumns));
                wei = find(cellfun(@(x) ~isempty(x),strfind(lower(app.SelectdataweightscolumnDropDown.Items),'weight')));
                if ~isempty(wei)
                    app.SelectdataweightscolumnDropDown.Value = app.SelectdataweightscolumnDropDown.Items{wei(1)};
                else
                    app.SelectdataweightscolumnDropDown.Value = app.SelectdataweightscolumnDropDown.Items{1};
                end
                
                app.dataWeightColumn = app.SelectdataweightscolumnDropDown.Value;
                
                app.SelectdatacolumnDropDown.Enable = true;
                app.SelectdatacolumnDropDownLabel.Enable = true;
                app.SelectdatacolumnDropDown.Items = sort(app.dataTable.Properties.VariableNames(numericColumns));
                dat = find(cellfun(@(x) ~isempty(x),strfind(lower(app.SelectdatacolumnDropDown.Items),'data')));
                if ~isempty(dat)
                    app.SelectdatacolumnDropDown.Value = app.SelectdatacolumnDropDown.Items{dat(1)};
                else
                    app.SelectdatacolumnDropDown.Value = app.SelectdatacolumnDropDown.Items{1};
                end
                
                app.dataColumn = app.SelectdatacolumnDropDown.Value;
                
                app.uselog10ofdataCheckBox.Enable = true;
                app.GeneratedatavectorsButton.Enable = true;
                
                app.UITable2.Data = []; drawnow;
                uialert(app.AnyPetroUIFigure,'Data loading finished.','Data Import',"Icon","success");
                app.CenterPanel.Enable = "off";
            end
            figure(app.AnyPetroUIFigure)
        end

        % Button pushed function: VieworeditdatatableButton
        function VieworeditdatatableButtonPushed(app, event)
            app.DialogApp = AnyPetro_displayDataTable(app,app.dataTable);
        end

        % Value changed function: 
        % SelectcolumnswithsamplespecificsListBox
        function SelectcolumnswithsamplespecificsListBoxValueChanged(app, event)
            app.sampleSpecifics = app.SelectcolumnswithsamplespecificsListBox.Value;
            for i = 1:numel(app.sampleSpecifics)
                %                 app.options.(app.sampleSpecifics{i}) = app.
            end
        end

        % Value changed function: SelectdataweightscolumnDropDown
        function SelectdataweightscolumnDropDownValueChanged(app, event)
            app.dataWeightColumn = app.SelectdataweightscolumnDropDown.Value;
        end

        % Close request function: AnyPetroUIFigure
        function AnyPetroUIFigureCloseRequest(app, event)
            delete(app.DialogApp)
            delete(app)
        end

        % Button pushed function: GeneratedatavectorsButton
        function GeneratedatavectorsButtonPushed(app, event)
            
            app.dobs = app.dataTable.(app.dataColumn);
            app.dobs = app.dobs(logical(app.dataTable.(app.SelectcolumnfordatauseDropDown.Value)));
            
            if app.uselog10ofdataCheckBox.Value
                app.dobs = log10(app.dobs);
            end
            
            app.err = app.dataTable.(app.dataWeightColumn);
            app.err = app.err(logical(app.dataTable.(app.SelectcolumnfordatauseDropDown.Value)));
            
            app.UITable.Data = [app.dobs app.err];
            app.UITable.RowName = 'numbered';
            app.UITable.Enable = 'on';
            
            if sum(isnan(app.UITable.Data),[1 2]) > 0
                app.CommandlineoutputTextArea.Value = 'Some data or data weights are missing (NaN). Inversion will not work.';
            else
                app.CommandlineoutputTextArea.Value = '';
            end
            
            app.options.DataType = app.dataTable.(app.SelectdatatypecolumnDropDown.Value);
            app.options.DataType = app.options.DataType(logical(app.dataTable.(app.SelectcolumnfordatauseDropDown.Value)));
            
            for i = 1:numel(app.sampleSpecifics)
                app.options.(app.sampleSpecifics{i}) = app.dataTable(logical(app.dataTable.(app.SelectcolumnfordatauseDropDown.Value)),:).(app.sampleSpecifics{i});
            end
            app.StartinversionButton.Enable = "off";
            app.UITable2.Data = []; drawnow;
            app.CenterPanel.Enable = "on";
            figure(app.AnyPetroUIFigure)
        end

        % Value changed function: SelectdatacolumnDropDown
        function SelectdatacolumnDropDownValueChanged(app, event)
            app.dataColumn = app.SelectdatacolumnDropDown.Value;
            figure(app.AnyPetroUIFigure)
        end

        % Button pushed function: SelectForwardOperatorButton
        function SelectForwardOperatorButtonPushed(app, event)
            if isempty(app.FOfilePath)
                if isempty(app.dataFilePath)
                    [FOfile,app.FOfilePath,~] = uigetfile({'*.txt'},'Select forward operator');
                else
                    [FOfile,app.FOfilePath,~] = uigetfile({'*.txt'},'Select forward operator',app.dataFilePath);
                end
            else
                [FOfile,app.FOfilePath,~] = uigetfile({'*.txt'},'Select forward operator',app.FOfilePath);
            end
            
            if ~isequal(FOfile,0)
                [ ~, app.FOname,~] = fileparts([app.FOfilePath,FOfile]);
                
                fid = fopen([app.FOfilePath FOfile]);
                T   = textscan(fid,'%s','delimiter','\n');
                fclose(fid);
                T = T{1,1};
                app.options.FOfile = T;
                
                idxName = find(strcmp(T,'[ModelNameStart]'))+1;
                idxMod = [find(strcmp(T,'[ModelParametersStart]'))+1; ...
                    find(strcmp(T,'[ModelParametersEnd]'))-1];
                opts = detectImportOptions([app.FOfilePath FOfile],'NumHeaderLines',idxMod(1));
                opts.DataLines = [idxMod(1)+1 idxMod(2)];
                opts.ConsecutiveDelimitersRule = 'join';
                opts.ExtraColumnsRule = 'ignore';
                opts.VariableNamesLine = idxMod(1);
                opts.Delimiter = '\t';
                opts.VariableNames = ["Name","lowerBound","upperBound","startingValue","referenceValue","weight","applyC1C2"];
                opts.SelectedVariableNames = ["Name","lowerBound","upperBound","startingValue","referenceValue","weight","applyC1C2"];
                opts.VariableTypes         = {'char', 'double', 'double', 'double', 'double', 'double', 'double'};
                app.modelTable = readtable([app.FOfilePath FOfile],opts);
                app.modelTable = app.modelTable(:,ismember(app.modelTable.Properties.VariableNames, ...
                    {'Name','lowerBound','upperBound','startingValue','referenceValue','weight','applyC1C2'}));
                
                for i = 2:size(app.modelTable,2)
                    try
                        app.modelTable.(app.modelTable.Properties.VariableNames{i}) = ...
                            str2double(app.modelTable.(app.modelTable.Properties.VariableNames{i}));
                    end
                end
                
                app.NoModPar = size(app.modelTable,1);
                app.FOLabel.Text = {T{idxName}; [num2str(app.NoModPar) ' model parameters']};
                app.FOLabel.Visible = 'on';
                app.ModelParameterNames = app.modelTable.Name;
                
                idxAux = [find(strcmp(T,'[AuxiliaryStatementsStart]'))+1; ...
                    find(strcmp(T,'[AuxiliaryStatementsEnd]'))-1];
                if ~isempty(idxAux)
                    app.options.AuxiliaryStatements = T(min(idxAux):max(idxAux));
                    % Execute Aux once, in case anonymous functions are
                    % included, that are used by the FO later. Those have to be
                    % defined when the FO is initialized.
                    for i = 1:numel(app.options.AuxiliaryStatements)
                        eval(app.options.AuxiliaryStatements{i});
                    end
                else
                    if isfield(app.options,'AuxiliaryStatements')
                        app.options = rmfield(app.options,'AuxiliaryStatements');
                    end
                end
                
                idxDat = unique([find(strcmp(T,'[SyntheticDataCalculationStart]'))+1; ...
                    find(strcmp(T,'[SyntheticDataCalculationEnd]'))-1]);
                opts = detectImportOptions([app.FOfilePath FOfile],'NumHeaderLines',idxDat(1));
                opts.DataLines = [idxDat(1)+1 idxDat(2)];
                opts.VariableNamesLine = idxDat(1);
                opts.Delimiter = '\t';
                opts.ConsecutiveDelimitersRule = 'join';
                opts.ExtraColumnsRule = 'ignore';
                opts.VariableNames = ["DataType", "Expression"];
                opts.SelectedVariableNames = ["DataType", "Expression"];
                opts.VariableTypes         = {'double', 'char'};
                
                dataExpr = readtable([app.FOfilePath FOfile],opts);
                
                idxSam = find(strcmp(T,'[SampleSpecificsStart]'))+1;
                if ~isempty(idxSam)
                    SamInput = T(idxSam);
                    SamInput = split(SamInput,',');
                    FOinputTemp = 'app.options.FOinput = @(options,i) [';
                    for i = 1:numel(SamInput)
                        FOinputTemp = [FOinputTemp 'options.' SamInput{i} '(i); '];
                    end
                    FOinputTemp = [FOinputTemp '];'];
                    eval(FOinputTemp);
                else
                    SamInput = {};
                    if isfield(app.options,'FOinput')
                        app.options = rmfield(app.options,'FOinput');
                    end
                end
                
                idxAdd = find(strcmp(T,'[AdditionalInputStart]'))+1;
                if ~isempty(idxAdd)
                    AddInput = T(idxAdd);
                    app.options.AddInput = AddInput;
                else
                    AddInput = {};
                    if isfield(app.options,'AddInput')
                        app.options = rmfield(app.options,'AddInput');
                    end
                end
                
                app.options.FO = {};
                for i=1:size(dataExpr,1)
                    if ~isempty(SamInput) && ~isempty(AddInput)
                        eval(cell2mat(['app.options.FO{i} = @(mod,' join([AddInput,SamInput],',') ') ' dataExpr(i,:).Expression ';']))
                    elseif ~isempty(SamInput) && isempty(AddInput)
                        eval(cell2mat(['app.options.FO{i} = @(mod,' join(SamInput,',') ') ' dataExpr(i,:).Expression ';']))
                    elseif isempty(SamInput) && ~isempty(AddInput)
                        eval(cell2mat(['app.options.FO{i} = @(mod,' join(AddInput,',') ') ' dataExpr(i,:).Expression ';']))
                    else
                        eval(cell2mat(['app.options.FO{i} = @(mod) ' dataExpr(i,:).Expression ';']))
                    end
                end
                
                app.UITable2.Data = app.modelTable(:,2:end);
                app.UITable2.RowName = app.ModelParameterNames;
                app.StartinversionButton.Enable = "on";
            end
            figure(app.AnyPetroUIFigure)
        end

        % Button pushed function: StartinversionButton
        function StartinversionButtonPushed(app, event)
            DoInversion(app)
        end

        % Button pushed function: ClearcommandlineButton
        function ClearcommandlineButtonPushed(app, event)
            app.CommandlineoutputTextArea.Value = '';
        end

        % Value changed function: CommandlineoutputTextArea
        function CommandlineoutputTextAreaValueChanged(app, event)
            app.CommandlineoutputTextArea.scroll('bottom');
            drawnow;
        end

        % Button pushed function: ExportresultButton
        function ExportresultButtonPushed(app, event)
            selpath = uigetdir;
            
            if ~isequal(selpath,0)
                exportDataTable = app.dataTable;
                exportDataTable = addvars(exportDataTable,nan(size(exportDataTable,1),1),'NewVariableNames','calculatedData');
                exportDataTable.calculatedData(exportDataTable.(app.SelectcolumnfordatauseDropDown.Value)) = app.data_calc_final;
                
                exportModelTable = app.modelTable;
                exportModelTable = addvars(exportModelTable,app.mod_erg,'NewVariableNames','inversionResult');
                exportModelTable = addvars(exportModelTable,app.mod_std,'NewVariableNames','parameterSTD');
                
                myTime = datestr(now,'yyyy-mm-dd_HHMMSS');
                writetable(exportDataTable, [selpath '\AnyPetro_inversionResult_' myTime '_data.txt'], "FileType", 'text','delimiter','\t');
                writetable(exportDataTable, [selpath '\AnyPetro_inversionResult_' myTime '_data.xlsx'], "FileType", 'spreadsheet');
                writetable(exportModelTable, [selpath '\AnyPetro_inversionResult_' myTime '_model.txt'], "FileType", 'text','delimiter','\t');
                writetable(exportModelTable, [selpath '\AnyPetro_inversionResult_' myTime '_model.xlsx'], "FileType", 'spreadsheet');
                
                exportgraphics(app.UIAxes,[selpath '\AnyPetro_inversionResult_' myTime '_Minimization.jpg'], 'Resolution', 300)
                exportgraphics(app.UIAxes2,[selpath '\AnyPetro_inversionResult_' myTime '_dataFit.jpg'], 'Resolution', 300)
                
                fid = fopen([selpath '\AnyPetro_inversionResult_' myTime '_log.txt'],'w');
                fprintf(fid,'%s\n','------------------');
                fprintf(fid,'%s\n','Forward Operator');
                fprintf(fid,'%s\n','------------------');
                fprintf(fid,'%s\n',app.options.FOfile{:});
                fprintf(fid,'%s\n','------------------');
                fprintf(fid,'%s\n','------------------');
                fprintf(fid,'%s\n',app.CommandlineoutputTextArea.Value{:});
                fclose(fid);
                uialert(app.AnyPetroUIFigure,'Export finished.','Export',"Icon","success");
            end
            figure(app.AnyPetroUIFigure)
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create AnyPetroUIFigure and hide until all components are created
            app.AnyPetroUIFigure = uifigure('Visible', 'off');
            app.AnyPetroUIFigure.AutoResizeChildren = 'off';
            app.AnyPetroUIFigure.Position = [100 100 1139 610];
            app.AnyPetroUIFigure.Name = 'AnyPetro';
            app.AnyPetroUIFigure.CloseRequestFcn = createCallbackFcn(app, @AnyPetroUIFigureCloseRequest, true);
            app.AnyPetroUIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.AnyPetroUIFigure);
            app.GridLayout.ColumnWidth = {362, '1x', 443};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.TitlePosition = 'centertop';
            app.LeftPanel.Title = 'Data Input';
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create GridLayout3
            app.GridLayout3 = uigridlayout(app.LeftPanel);
            app.GridLayout3.RowHeight = {'1x', '1x', '1x', '1x', '1x', '1x', '3x', '2x', '2x'};

            % Create LoaddatafileButton
            app.LoaddatafileButton = uibutton(app.GridLayout3, 'push');
            app.LoaddatafileButton.ButtonPushedFcn = createCallbackFcn(app, @LoaddatafileButtonPushed, true);
            app.LoaddatafileButton.Layout.Row = 1;
            app.LoaddatafileButton.Layout.Column = 1;
            app.LoaddatafileButton.Text = 'Load data file';

            % Create VieworeditdatatableButton
            app.VieworeditdatatableButton = uibutton(app.GridLayout3, 'push');
            app.VieworeditdatatableButton.ButtonPushedFcn = createCallbackFcn(app, @VieworeditdatatableButtonPushed, true);
            app.VieworeditdatatableButton.Enable = 'off';
            app.VieworeditdatatableButton.Layout.Row = 1;
            app.VieworeditdatatableButton.Layout.Column = 2;
            app.VieworeditdatatableButton.Text = 'View or edit data table';

            % Create GeneratedatavectorsButton
            app.GeneratedatavectorsButton = uibutton(app.GridLayout3, 'push');
            app.GeneratedatavectorsButton.ButtonPushedFcn = createCallbackFcn(app, @GeneratedatavectorsButtonPushed, true);
            app.GeneratedatavectorsButton.Enable = 'off';
            app.GeneratedatavectorsButton.Layout.Row = 6;
            app.GeneratedatavectorsButton.Layout.Column = 1;
            app.GeneratedatavectorsButton.Text = 'Generate data vectors';

            % Create uselog10ofdataCheckBox
            app.uselog10ofdataCheckBox = uicheckbox(app.GridLayout3);
            app.uselog10ofdataCheckBox.Enable = 'off';
            app.uselog10ofdataCheckBox.Text = 'use log10 of data';
            app.uselog10ofdataCheckBox.Layout.Row = 6;
            app.uselog10ofdataCheckBox.Layout.Column = 2;
            app.uselog10ofdataCheckBox.Value = true;

            % Create SelectdataweightscolumnDropDown
            app.SelectdataweightscolumnDropDown = uidropdown(app.GridLayout3);
            app.SelectdataweightscolumnDropDown.Items = {};
            app.SelectdataweightscolumnDropDown.ValueChangedFcn = createCallbackFcn(app, @SelectdataweightscolumnDropDownValueChanged, true);
            app.SelectdataweightscolumnDropDown.Enable = 'off';
            app.SelectdataweightscolumnDropDown.Layout.Row = 5;
            app.SelectdataweightscolumnDropDown.Layout.Column = 2;
            app.SelectdataweightscolumnDropDown.Value = {};

            % Create SelectdataweightscolumnDropDownLabel
            app.SelectdataweightscolumnDropDownLabel = uilabel(app.GridLayout3);
            app.SelectdataweightscolumnDropDownLabel.HorizontalAlignment = 'right';
            app.SelectdataweightscolumnDropDownLabel.Enable = 'off';
            app.SelectdataweightscolumnDropDownLabel.Layout.Row = 5;
            app.SelectdataweightscolumnDropDownLabel.Layout.Column = 1;
            app.SelectdataweightscolumnDropDownLabel.Text = 'Select data weights column';

            % Create SelectcolumnswithsamplespecificsListBox
            app.SelectcolumnswithsamplespecificsListBox = uilistbox(app.GridLayout3);
            app.SelectcolumnswithsamplespecificsListBox.Items = {};
            app.SelectcolumnswithsamplespecificsListBox.Multiselect = 'on';
            app.SelectcolumnswithsamplespecificsListBox.ValueChangedFcn = createCallbackFcn(app, @SelectcolumnswithsamplespecificsListBoxValueChanged, true);
            app.SelectcolumnswithsamplespecificsListBox.Enable = 'off';
            app.SelectcolumnswithsamplespecificsListBox.Layout.Row = 7;
            app.SelectcolumnswithsamplespecificsListBox.Layout.Column = 2;
            app.SelectcolumnswithsamplespecificsListBox.Value = {};

            % Create SelectcolumnswithsamplespecificsListBoxLabel
            app.SelectcolumnswithsamplespecificsListBoxLabel = uilabel(app.GridLayout3);
            app.SelectcolumnswithsamplespecificsListBoxLabel.HorizontalAlignment = 'right';
            app.SelectcolumnswithsamplespecificsListBoxLabel.WordWrap = 'on';
            app.SelectcolumnswithsamplespecificsListBoxLabel.Enable = 'off';
            app.SelectcolumnswithsamplespecificsListBoxLabel.Layout.Row = 7;
            app.SelectcolumnswithsamplespecificsListBoxLabel.Layout.Column = 1;
            app.SelectcolumnswithsamplespecificsListBoxLabel.Text = 'Select columns with sample specifics';

            % Create UITable
            app.UITable = uitable(app.GridLayout3);
            app.UITable.ColumnName = {'data'; 'weight'};
            app.UITable.RowName = {};
            app.UITable.Enable = 'off';
            app.UITable.Layout.Row = [8 9];
            app.UITable.Layout.Column = [1 2];

            % Create SelectdatatypecolumnDropDownLabel
            app.SelectdatatypecolumnDropDownLabel = uilabel(app.GridLayout3);
            app.SelectdatatypecolumnDropDownLabel.HorizontalAlignment = 'right';
            app.SelectdatatypecolumnDropDownLabel.WordWrap = 'on';
            app.SelectdatatypecolumnDropDownLabel.Enable = 'off';
            app.SelectdatatypecolumnDropDownLabel.Layout.Row = 3;
            app.SelectdatatypecolumnDropDownLabel.Layout.Column = 1;
            app.SelectdatatypecolumnDropDownLabel.Text = 'Select data type column';

            % Create SelectdatatypecolumnDropDown
            app.SelectdatatypecolumnDropDown = uidropdown(app.GridLayout3);
            app.SelectdatatypecolumnDropDown.Items = {};
            app.SelectdatatypecolumnDropDown.Enable = 'off';
            app.SelectdatatypecolumnDropDown.Layout.Row = 3;
            app.SelectdatatypecolumnDropDown.Layout.Column = 2;
            app.SelectdatatypecolumnDropDown.Value = {};

            % Create SelectcolumnfordatauseDropDown
            app.SelectcolumnfordatauseDropDown = uidropdown(app.GridLayout3);
            app.SelectcolumnfordatauseDropDown.Items = {};
            app.SelectcolumnfordatauseDropDown.Enable = 'off';
            app.SelectcolumnfordatauseDropDown.Layout.Row = 2;
            app.SelectcolumnfordatauseDropDown.Layout.Column = 2;
            app.SelectcolumnfordatauseDropDown.Value = {};

            % Create SelectcolumnfordatauseDropDownLabel
            app.SelectcolumnfordatauseDropDownLabel = uilabel(app.GridLayout3);
            app.SelectcolumnfordatauseDropDownLabel.HorizontalAlignment = 'right';
            app.SelectcolumnfordatauseDropDownLabel.WordWrap = 'on';
            app.SelectcolumnfordatauseDropDownLabel.Enable = 'off';
            app.SelectcolumnfordatauseDropDownLabel.Layout.Row = 2;
            app.SelectcolumnfordatauseDropDownLabel.Layout.Column = 1;
            app.SelectcolumnfordatauseDropDownLabel.Text = 'Select column for data use';

            % Create SelectdatacolumnDropDownLabel
            app.SelectdatacolumnDropDownLabel = uilabel(app.GridLayout3);
            app.SelectdatacolumnDropDownLabel.HorizontalAlignment = 'right';
            app.SelectdatacolumnDropDownLabel.Enable = 'off';
            app.SelectdatacolumnDropDownLabel.Layout.Row = 4;
            app.SelectdatacolumnDropDownLabel.Layout.Column = 1;
            app.SelectdatacolumnDropDownLabel.Text = 'Select data column';

            % Create SelectdatacolumnDropDown
            app.SelectdatacolumnDropDown = uidropdown(app.GridLayout3);
            app.SelectdatacolumnDropDown.Items = {};
            app.SelectdatacolumnDropDown.ValueChangedFcn = createCallbackFcn(app, @SelectdatacolumnDropDownValueChanged, true);
            app.SelectdatacolumnDropDown.Enable = 'off';
            app.SelectdatacolumnDropDown.Layout.Row = 4;
            app.SelectdatacolumnDropDown.Layout.Column = 2;
            app.SelectdatacolumnDropDown.Value = {};

            % Create CenterPanel
            app.CenterPanel = uipanel(app.GridLayout);
            app.CenterPanel.Enable = 'off';
            app.CenterPanel.TitlePosition = 'centertop';
            app.CenterPanel.Title = 'Inversion Set-up';
            app.CenterPanel.Layout.Row = 1;
            app.CenterPanel.Layout.Column = 2;

            % Create GridLayout2
            app.GridLayout2 = uigridlayout(app.CenterPanel);
            app.GridLayout2.RowHeight = {'1x', '3x', '2x', '1x', '1x'};

            % Create UITable2
            app.UITable2 = uitable(app.GridLayout2);
            app.UITable2.ColumnName = {'lower bound'; 'upper bound'; 'starting value'; 'reference value'; 'weight'; 'applyC1C2'};
            app.UITable2.RowName = {};
            app.UITable2.ColumnEditable = true;
            app.UITable2.Layout.Row = 2;
            app.UITable2.Layout.Column = [1 2];

            % Create SelectForwardOperatorButton
            app.SelectForwardOperatorButton = uibutton(app.GridLayout2, 'push');
            app.SelectForwardOperatorButton.ButtonPushedFcn = createCallbackFcn(app, @SelectForwardOperatorButtonPushed, true);
            app.SelectForwardOperatorButton.Layout.Row = 1;
            app.SelectForwardOperatorButton.Layout.Column = 1;
            app.SelectForwardOperatorButton.Text = 'Select Forward Operator';

            % Create FOLabel
            app.FOLabel = uilabel(app.GridLayout2);
            app.FOLabel.HorizontalAlignment = 'center';
            app.FOLabel.WordWrap = 'on';
            app.FOLabel.Visible = 'off';
            app.FOLabel.Layout.Row = 1;
            app.FOLabel.Layout.Column = 2;
            app.FOLabel.Text = 'FO';

            % Create RegularizationPanel
            app.RegularizationPanel = uipanel(app.GridLayout2);
            app.RegularizationPanel.TitlePosition = 'centertop';
            app.RegularizationPanel.Title = 'Regularization';
            app.RegularizationPanel.Layout.Row = 3;
            app.RegularizationPanel.Layout.Column = 1;

            % Create GridLayout5
            app.GridLayout5 = uigridlayout(app.RegularizationPanel);
            app.GridLayout5.ColumnWidth = {'1x', '2x'};
            app.GridLayout5.RowHeight = {'1x', '1x', '1x'};

            % Create lambda0EditFieldLabel
            app.lambda0EditFieldLabel = uilabel(app.GridLayout5);
            app.lambda0EditFieldLabel.HorizontalAlignment = 'right';
            app.lambda0EditFieldLabel.Layout.Row = 1;
            app.lambda0EditFieldLabel.Layout.Column = 1;
            app.lambda0EditFieldLabel.Text = 'lambda0';

            % Create lambda0EditField
            app.lambda0EditField = uieditfield(app.GridLayout5, 'numeric');
            app.lambda0EditField.Limits = [0 Inf];
            app.lambda0EditField.ValueDisplayFormat = '%.2e';
            app.lambda0EditField.Layout.Row = 1;
            app.lambda0EditField.Layout.Column = 2;

            % Create lambda1EditFieldLabel
            app.lambda1EditFieldLabel = uilabel(app.GridLayout5);
            app.lambda1EditFieldLabel.HorizontalAlignment = 'right';
            app.lambda1EditFieldLabel.Layout.Row = 2;
            app.lambda1EditFieldLabel.Layout.Column = 1;
            app.lambda1EditFieldLabel.Text = 'lambda1';

            % Create lambda1EditField
            app.lambda1EditField = uieditfield(app.GridLayout5, 'numeric');
            app.lambda1EditField.Limits = [0 Inf];
            app.lambda1EditField.ValueDisplayFormat = '%.2e';
            app.lambda1EditField.Layout.Row = 2;
            app.lambda1EditField.Layout.Column = 2;

            % Create lambda2EditFieldLabel
            app.lambda2EditFieldLabel = uilabel(app.GridLayout5);
            app.lambda2EditFieldLabel.HorizontalAlignment = 'right';
            app.lambda2EditFieldLabel.Layout.Row = 3;
            app.lambda2EditFieldLabel.Layout.Column = 1;
            app.lambda2EditFieldLabel.Text = 'lambda2';

            % Create lambda2EditField
            app.lambda2EditField = uieditfield(app.GridLayout5, 'numeric');
            app.lambda2EditField.Limits = [0 Inf];
            app.lambda2EditField.ValueDisplayFormat = '%.2e';
            app.lambda2EditField.Layout.Row = 3;
            app.lambda2EditField.Layout.Column = 2;

            % Create InversionParametersPanel
            app.InversionParametersPanel = uipanel(app.GridLayout2);
            app.InversionParametersPanel.TitlePosition = 'centertop';
            app.InversionParametersPanel.Title = 'Inversion Parameters';
            app.InversionParametersPanel.Layout.Row = [4 5];
            app.InversionParametersPanel.Layout.Column = 1;

            % Create GridLayout6
            app.GridLayout6 = uigridlayout(app.InversionParametersPanel);
            app.GridLayout6.RowHeight = {'1x', '1x', '1x'};

            % Create MaximumnumberofiterationsEditFieldLabel
            app.MaximumnumberofiterationsEditFieldLabel = uilabel(app.GridLayout6);
            app.MaximumnumberofiterationsEditFieldLabel.HorizontalAlignment = 'right';
            app.MaximumnumberofiterationsEditFieldLabel.WordWrap = 'on';
            app.MaximumnumberofiterationsEditFieldLabel.Layout.Row = 1;
            app.MaximumnumberofiterationsEditFieldLabel.Layout.Column = 1;
            app.MaximumnumberofiterationsEditFieldLabel.Text = 'Maximum number of iterations';

            % Create MaximumnumberofiterationsEditField
            app.MaximumnumberofiterationsEditField = uieditfield(app.GridLayout6, 'numeric');
            app.MaximumnumberofiterationsEditField.Limits = [1 Inf];
            app.MaximumnumberofiterationsEditField.RoundFractionalValues = 'on';
            app.MaximumnumberofiterationsEditField.ValueDisplayFormat = '%.0f';
            app.MaximumnumberofiterationsEditField.Layout.Row = 1;
            app.MaximumnumberofiterationsEditField.Layout.Column = 2;
            app.MaximumnumberofiterationsEditField.Value = 10;

            % Create PerturbationEditFieldLabel
            app.PerturbationEditFieldLabel = uilabel(app.GridLayout6);
            app.PerturbationEditFieldLabel.HorizontalAlignment = 'right';
            app.PerturbationEditFieldLabel.WordWrap = 'on';
            app.PerturbationEditFieldLabel.Layout.Row = 2;
            app.PerturbationEditFieldLabel.Layout.Column = 1;
            app.PerturbationEditFieldLabel.Text = 'Perturbation';

            % Create PerturbationEditField
            app.PerturbationEditField = uieditfield(app.GridLayout6, 'numeric');
            app.PerturbationEditField.Limits = [1e-20 Inf];
            app.PerturbationEditField.ValueDisplayFormat = '%.2e';
            app.PerturbationEditField.Layout.Row = 2;
            app.PerturbationEditField.Layout.Column = 2;
            app.PerturbationEditField.Value = 0.0001;

            % Create ToleranceEditFieldLabel
            app.ToleranceEditFieldLabel = uilabel(app.GridLayout6);
            app.ToleranceEditFieldLabel.HorizontalAlignment = 'right';
            app.ToleranceEditFieldLabel.WordWrap = 'on';
            app.ToleranceEditFieldLabel.Layout.Row = 3;
            app.ToleranceEditFieldLabel.Layout.Column = 1;
            app.ToleranceEditFieldLabel.Text = 'Tolerance';

            % Create ToleranceEditField
            app.ToleranceEditField = uieditfield(app.GridLayout6, 'numeric');
            app.ToleranceEditField.Limits = [1e-20 Inf];
            app.ToleranceEditField.ValueDisplayFormat = '%.2e';
            app.ToleranceEditField.Layout.Row = 3;
            app.ToleranceEditField.Layout.Column = 2;
            app.ToleranceEditField.Value = 1e-05;

            % Create StartinversionButton
            app.StartinversionButton = uibutton(app.GridLayout2, 'push');
            app.StartinversionButton.ButtonPushedFcn = createCallbackFcn(app, @StartinversionButtonPushed, true);
            app.StartinversionButton.FontWeight = 'bold';
            app.StartinversionButton.Layout.Row = 4;
            app.StartinversionButton.Layout.Column = 2;
            app.StartinversionButton.Text = 'Start inversion';

            % Create ParametertransformationButtonGroup
            app.ParametertransformationButtonGroup = uibuttongroup(app.GridLayout2);
            app.ParametertransformationButtonGroup.TitlePosition = 'centertop';
            app.ParametertransformationButtonGroup.Title = 'Parameter transformation';
            app.ParametertransformationButtonGroup.Layout.Row = 3;
            app.ParametertransformationButtonGroup.Layout.Column = 2;

            % Create NoneButton
            app.NoneButton = uiradiobutton(app.ParametertransformationButtonGroup);
            app.NoneButton.Text = 'None';
            app.NoneButton.Position = [11 10 58 22];

            % Create LogButton
            app.LogButton = uiradiobutton(app.ParametertransformationButtonGroup);
            app.LogButton.Text = 'Log';
            app.LogButton.Position = [11 45 65 22];

            % Create RangeButton
            app.RangeButton = uiradiobutton(app.ParametertransformationButtonGroup);
            app.RangeButton.Text = 'Range';
            app.RangeButton.WordWrap = 'on';
            app.RangeButton.Position = [11 80 85 22];
            app.RangeButton.Value = true;

            % Create ExportresultButton
            app.ExportresultButton = uibutton(app.GridLayout2, 'push');
            app.ExportresultButton.ButtonPushedFcn = createCallbackFcn(app, @ExportresultButtonPushed, true);
            app.ExportresultButton.Enable = 'off';
            app.ExportresultButton.Layout.Row = 5;
            app.ExportresultButton.Layout.Column = 2;
            app.ExportresultButton.Text = 'Export result';

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.TitlePosition = 'centertop';
            app.RightPanel.Title = 'Output';
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 3;

            % Create GridLayout4
            app.GridLayout4 = uigridlayout(app.RightPanel);
            app.GridLayout4.RowHeight = {'1x', '5x', '6x'};

            % Create UIAxes
            app.UIAxes = uiaxes(app.GridLayout4);
            title(app.UIAxes, 'Minimization')
            xlabel(app.UIAxes, 'iteration')
            ylabel(app.UIAxes, 'normalized objective function')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.PlotBoxAspectRatio = [1 1 1];
            app.UIAxes.TickLabelInterpreter = 'latex';
            app.UIAxes.YScale = 'log';
            app.UIAxes.YMinorTick = 'on';
            app.UIAxes.XGrid = 'on';
            app.UIAxes.YGrid = 'on';
            app.UIAxes.YMinorGrid = 'on';
            app.UIAxes.Layout.Row = 3;
            app.UIAxes.Layout.Column = 1;

            % Create UIAxes2
            app.UIAxes2 = uiaxes(app.GridLayout4);
            title(app.UIAxes2, 'Data fit')
            xlabel(app.UIAxes2, 'measured')
            ylabel(app.UIAxes2, 'inversion')
            zlabel(app.UIAxes2, 'Z')
            app.UIAxes2.PlotBoxAspectRatio = [1 1 1];
            app.UIAxes2.TickLabelInterpreter = 'latex';
            app.UIAxes2.XGrid = 'on';
            app.UIAxes2.YGrid = 'on';
            app.UIAxes2.Layout.Row = 3;
            app.UIAxes2.Layout.Column = 2;

            % Create CommandlineoutputTextArea
            app.CommandlineoutputTextArea = uitextarea(app.GridLayout4);
            app.CommandlineoutputTextArea.ValueChangedFcn = createCallbackFcn(app, @CommandlineoutputTextAreaValueChanged, true);
            app.CommandlineoutputTextArea.Editable = 'off';
            app.CommandlineoutputTextArea.Layout.Row = 2;
            app.CommandlineoutputTextArea.Layout.Column = [1 2];

            % Create CommandlineoutputTextAreaLabel
            app.CommandlineoutputTextAreaLabel = uilabel(app.GridLayout4);
            app.CommandlineoutputTextAreaLabel.HorizontalAlignment = 'center';
            app.CommandlineoutputTextAreaLabel.Layout.Row = 1;
            app.CommandlineoutputTextAreaLabel.Layout.Column = 1;
            app.CommandlineoutputTextAreaLabel.Text = 'Command line output';

            % Create ClearcommandlineButton
            app.ClearcommandlineButton = uibutton(app.GridLayout4, 'push');
            app.ClearcommandlineButton.ButtonPushedFcn = createCallbackFcn(app, @ClearcommandlineButtonPushed, true);
            app.ClearcommandlineButton.Layout.Row = 1;
            app.ClearcommandlineButton.Layout.Column = 2;
            app.ClearcommandlineButton.Text = 'Clear command line';

            % Show the figure after all components are created
            app.AnyPetroUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = AnyPetro_reflow_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.AnyPetroUIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.AnyPetroUIFigure)
        end
    end
end