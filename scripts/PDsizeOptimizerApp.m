function PDsizeOptimizerApp
%PDSIZEOPTIMIZERAPP Displays data from retroreflector expermients
%   Run this app and it will plot pre-computed data dynamically, as well as
%   the theoretical results
    
%% init simulation parameters
D_d = .004:.002:.016;

%%   create GUI
    fig = uifigure;
    fig.Name = "Photodiode Optimizer";
    fig.WindowState = 'maximized';

%   create master grid layout
    gl = uigridlayout(fig,[4 3]);
    gl.RowHeight =  {50, '1x', '1x', 50};
    gl.ColumnWidth = {'1x', '2x', '2x'};
    
%   label
    dir_lbl = uilabel(gl);
    dir_lbl.Text = "Directions here blah blah blah blah blah";
    dir_lbl.Layout.Row = 1;
    dir_lbl.Layout.Column = 1;
    dir_lbl.WordWrap = 'on';

    title = uilabel(gl);
    title.Text = "Title";
    title.Layout.Row = 1;
    title.Layout.Column = [2 3];

    meta = uilabel(gl);
    meta.Text = "Meta Data";
    meta.Layout.Row = 4;
    meta.Layout.Column = [2 3];


%     Add tab group on left
    tabgp = uitabgroup(gl);
    tabgp.Layout.Row = [2 4];
    tabgp.Layout.Column = 1;
    
%%   Add tab for basic controls, fill in with widgets
    ctrl_tab = uitab(tabgp,'Title', 'Control');
    ctrl_tab_gl = uigridlayout(ctrl_tab,[10,3]);
    ctrl_tab_gl.ColumnWidth = {'1x', '1x', '1x'};
    ctrl_tab_gl.RowHeight = {'1x', 'fit'};
    
%     Select Diode Size
    options = split(sprintf('%dmm ', D_d*1000 ));
    options(end) = []; %remove null char

    pd_box = uilistbox(ctrl_tab_gl, ...
        'Items', options,...
        'ItemsData', D_d,...
        'Multiselect','on');
    pd_box.Layout.Row = [2 5];
    pd_box.Layout.Column = [1 3];

%%   Add advanced tab
    adv_tab = uitab(tabgp,'Title', 'Advanced');

%%    init all graphs
%   for alpha as function of horiz displacement
    ax_alpha_v = uiaxes(gl);
    ax_alpha_v.Layout.Row = 2;
    ax_alpha_v.Layout.Column = 2;

%   alpha as function of PD size
    ax_alpha_D_d = uiaxes(gl);
    ax_alpha_D_d.Layout.Row = 2;
    ax_alpha_D_d.Layout.Column = 3;

%   normalized power as function of horiz disp
    ax_pwr_v = uiaxes(gl);
    ax_pwr_v.Layout.Row = 3;
    ax_pwr_v.Layout.Column = 2;

%   normalized power as function of PD size
    ax_pwr_D_d = uiaxes(gl);
    ax_pwr_D_d.Layout.Row = 3;
    ax_pwr_D_d.Layout.Column = 3;

%   plot defaults
    surf(ax_alpha_v,peaks);
    mesh(ax_alpha_D_d,peaks);
    waterfall(ax_pwr_v,peaks);
    surf(ax_pwr_D_d,peaks);
end


function changePlotType(src,event,ax)
    type = event.Value;
    switch type
        case "Surf"
            surf(ax,peaks);
        case "Mesh"
            mesh(ax,peaks);
        case "Waterfall"
            waterfall(ax,peaks);
    end
end

