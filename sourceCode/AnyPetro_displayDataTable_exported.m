classdef AnyPetro_displayDataTable_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        AnyPetrodisplaydatatableUIFigure  matlab.ui.Figure
        GridLayout          matlab.ui.container.GridLayout
        SaveandcloseButton  matlab.ui.control.Button
        UITable             matlab.ui.control.Table
    end

    
    properties (Access = private)
        CallingApp % Description
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app, phyany_reflow, dataTable)
            app.CallingApp = phyany_reflow;
            app.UITable.ColumnName = dataTable.Properties.VariableNames;
            app.UITable.RowName = 'numbered';
            app.UITable.Data = dataTable;
        end

        % Close request function: AnyPetrodisplaydatatableUIFigure
        function AnyPetrodisplaydatatableUIFigureCloseRequest(app, event)
            updateDataTable(app.CallingApp,app.UITable.Data)
            delete(app)
        end

        % Button pushed function: SaveandcloseButton
        function SaveandcloseButtonPushed(app, event)
            updateDataTable(app.CallingApp,app.UITable.Data)
            delete(app)
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create AnyPetrodisplaydatatableUIFigure and hide until all components are created
            app.AnyPetrodisplaydatatableUIFigure = uifigure('Visible', 'off');
            app.AnyPetrodisplaydatatableUIFigure.Position = [100 100 640 480];
            app.AnyPetrodisplaydatatableUIFigure.Name = 'AnyPetro - display data table';
            app.AnyPetrodisplaydatatableUIFigure.CloseRequestFcn = createCallbackFcn(app, @AnyPetrodisplaydatatableUIFigureCloseRequest, true);
            app.AnyPetrodisplaydatatableUIFigure.Scrollable = 'on';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.AnyPetrodisplaydatatableUIFigure);
            app.GridLayout.ColumnWidth = {'1x', '3x'};
            app.GridLayout.RowHeight = {'1x', '7.5x'};

            % Create UITable
            app.UITable = uitable(app.GridLayout);
            app.UITable.ColumnName = {'Column 1'; 'Column 2'; 'Column 3'; 'Column 4'};
            app.UITable.RowName = {};
            app.UITable.ColumnSortable = true;
            app.UITable.ColumnEditable = true;
            app.UITable.Layout.Row = 2;
            app.UITable.Layout.Column = [1 2];

            % Create SaveandcloseButton
            app.SaveandcloseButton = uibutton(app.GridLayout, 'push');
            app.SaveandcloseButton.ButtonPushedFcn = createCallbackFcn(app, @SaveandcloseButtonPushed, true);
            app.SaveandcloseButton.WordWrap = 'on';
            app.SaveandcloseButton.Layout.Row = 1;
            app.SaveandcloseButton.Layout.Column = 1;
            app.SaveandcloseButton.Text = 'Save and close';

            % Show the figure after all components are created
            app.AnyPetrodisplaydatatableUIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = AnyPetro_displayDataTable_exported(varargin)

            runningApp = getRunningApp(app);

            % Check for running singleton app
            if isempty(runningApp)

                % Create UIFigure and components
                createComponents(app)

                % Register the app with App Designer
                registerApp(app, app.AnyPetrodisplaydatatableUIFigure)

                % Execute the startup function
                runStartupFcn(app, @(app)startupFcn(app, varargin{:}))
            else

                % Focus the running singleton app
                figure(runningApp.AnyPetrodisplaydatatableUIFigure)

                app = runningApp;
            end

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.AnyPetrodisplaydatatableUIFigure)
        end
    end
end