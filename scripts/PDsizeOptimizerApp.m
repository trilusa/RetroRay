function PDsizeOptimizerApp
%PDSIZEOPTIMIZERAPP Displays data from retroreflector expermients
%   Run this app and it will plot pre-computed data dynamically, as well as
%   the theoretical results
    
%% load simulation parameters/data
    data = load('monte-carlo-alpha_N200000_PDALL_L50mm(1).mat');

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
    ctrl_tab_gl.RowHeight = {'fit', 'fit'};
    
%     Select and display Diode Size
    pd_txt = uilabel(ctrl_tab_gl);
    options = split(sprintf('%dmm ', data.D_d*1000 ));
    options(end) = []; %remove null char
    pd_box = uilistbox(ctrl_tab_gl, ...
        'Items', options,...
        'ItemsData', 1:length(data.D_d),...
        'Multiselect','on',...
        'ValueChangedFcn', {@PDselectionChanged, pd_txt, data});
    pd_box.Layout.Row = [3 5];
    pd_box.Layout.Column = [1 3];
    pd_box.Value = pd_box.ItemsData;

    pd_txt.Text = ['Selected Diode Sizes (mm): ' mat2str(data.D_d(pd_box.Value))];
    pd_txt.Layout.Column = [1 3];

%     Update button
    btn = uibutton(ctrl_tab_gl, 'push',...
        'Text', 'Plot');
%          'ButtonPushedFcn', {@UpdatePlots, data, pd_box, ax_alpha_v});
    btn.Layout.Row = 6;
    btn.Layout.Column = [1 3];


%%    init all graphs
%   for alpha as function of horiz displacement
    ax_alpha_v = uiaxes(gl);
    ax_alpha_v.Layout.Row = 2;
    ax_alpha_v.Layout.Column = 2;
    ax_alpha_v.Title.String = "Alpha as function of horiz displacement";
%     ax_alpha_v.Title = "Alpha as function of horiz displacement";
    ax_alpha_v.XLabel.String = 'v (m)';
    ax_alpha_v.YLabel.String = '\alpha';
    ax_alpha_v.XGrid = 'on';
    ax_alpha_v.YGrid = 'on';


%   alpha as function of PD size
    ax_alpha_D_d = uiaxes(gl);
    ax_alpha_D_d.Layout.Row = 2;
    ax_alpha_D_d.Layout.Column = 3;
    ax_alpha_D_d.Title.String = "alpha as function of PD size";
    ax_alpha_D_d.XLabel.String = 'PD diameter (mm)';
    ax_alpha_D_d.YLabel.String = 'Normlized Power';
    ax_alpha_D_d.XGrid = 'on';
    ax_alpha_D_d.YGrid = 'on';
    ax_alpha_D_d.XLim = [0 20];

%   normalized power as function of horiz disp
    ax_pwr_v = uiaxes(gl);
    ax_pwr_v.Layout.Row = 3;
    ax_pwr_v.Layout.Column = 2;
%     myplot(ax_alpha_v,  pd_box.Value, data.v, data.Pr_norm, data.data_cos_d_norm );
    ax_pwr_v.Title.String = "normalized power as function of horiz disp";
    ax_pwr_v.XLabel.String = 'v (m)';
    ax_pwr_v.YLabel.String = 'Normlized Power';
    ax_pwr_v.XGrid = 'on';
    ax_pwr_v.YGrid = 'on';
    

%   normalized power as function of PD size
    ax_pwr_D_d = uiaxes(gl);
    ax_pwr_D_d.Layout.Row = 3;
    ax_pwr_D_d.Layout.Column = 3;
    ax_pwr_D_d.Title.String = "normalized power as function of PD size";
    ax_pwr_D_d.XLabel.String = 'PD diameter (mm)';
    ax_pwr_D_d.YLabel.String = 'Normlized Power';
    ax_pwr_D_d.XGrid = 'on';
    ax_pwr_D_d.YGrid = 'on';
    ax_pwr_D_d.XLim = [0 20];

    ax = [ax_alpha_v, ax_alpha_D_d, ax_pwr_v, ax_pwr_D_d];
    btn.ButtonPushedFcn = {@UpdatePlots, data, pd_box, ax};
end

function PDselectionChanged(src,event,PDtxt,data)
    PDtxt.Text = ['Selected Diode Sizes (mm):' mat2str(data.D_d(event.Value))];
end

function UpdatePlots(src,event,data,pd_box,ax)
    src.Text = 'Update';
    keys = sort(pd_box.Value);
    D_d = data.D_d(keys);
    v = data.v;
    pwr_theory = data.Pr_norm(keys,:);
    pwr_sim = data.data_cos_d_norm(keys,:);
    pwr_error = (pwr_theory - pwr_sim) ./ pwr_theory;  

    plot(ax(3), v, [pwr_theory; pwr_sim ; pwr_error] );
    legend(ax(3),[num2str(D_d'*1000)  repmat(' mm',length(D_d),1)]);

    plot(ax(4), D_d*1000, [data.Pr_norm(keys,:)'; data.data_cos_d_norm(keys,:)'], 'o');
    legend(ax(4),[num2str(v')  repmat(' m',length(v),1)]);
   
end

function myplot(ax, keys, labels, x, theory, sim)

       plot(ax, x,[theory; sim]);

           %[data.data_cos_d_norm(pd_box.Value,:); data.Pr_norm(pd_box.Value,:)],['o', '-']);

%     plot(ax, data.v, data.data_cos_d_norm(val(1),:), 'x', 'SeriesIndex', 1);
%     hold(ax,'on')
%     plot(ax, data.v, data.Pr_norm(val(1),:), '-', 'SeriesIndex', 1);
%         
%     for i=2:length(val)
%         plot(ax, data.v, data.data_cos_d_norm(val(i),:), 'x', 'SeriesIndex', i);
%         plot(ax, data.v, data.Pr_norm(val(i),:), '-', 'SeriesIndex', i);
%     end
%     legend
%     hold(ax_alpha_v,'off')
end