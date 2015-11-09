function Happy_Birthday
handles.GUI_f=figure('visible','on','Position',[10,50,1350,630],'color',[1 1 1]);
h = uitabgroup('TabLocation','Top','Margin',2); drawnow;
t1=uitab(h,'title','Hello');
t2=uitab(h, 'title', 'My Cake');
t3=uitab(h, 'title', 'Bye');
%% Buttons for tab number 1
handles.title1=uicontrol(t1,'Style','edit','string','In the name of God',...
    'Position',[150,540,1000,50],'BackgroundColor',[1 1 1],'FontSize',18,'FontName',...
    'Eras Demi ITC');
A={'Hello Nooshin. First of all happy birthday (and also happy new year).'...
 'I wish you have a very good year and be successful anytime, anywhere.',...
'Nooshin, I''ve made a cake for you. I think it''s delicious and I hope you like it too.',...
'Go to page "My Cake". When you are ready click "Show my cake" button.'...
'The cake will be shown then and the "Candles" & "Blow off" buttons are appeared. '...
'Now click on "Candles" button to turn on the candles.'...
'First wish, then press on "Blow off" button to blow off the candles.'...
'Enjoy your cake & once again . . .'...
' "HAPPY BIRTHDAY" !!!'};
handles.title2=uicontrol(t1,'Style','edit','string',A(1)...
    ,'Position',[150,490,1000,50],'BackgroundColor',[1 1 1],'FontSize',15,'FontName',...
    'Eras Demi ITC');
handles.title3=uicontrol(t1,'Style','edit','string',A(2)...
    ,'Position',[150,440,1000,50],'BackgroundColor',[1 1 1],'FontSize',15,'FontName',...
    'Eras Demi ITC');
handles.title4=uicontrol(t1,'Style','edit','string',A(3)...
    ,'Position',[150,390,1000,50],'BackgroundColor',[1 1 1],'FontSize',15,'FontName',...
    'Eras Demi ITC');
handles.title5=uicontrol(t1,'Style','edit','string',A(4)...
    ,'Position',[150,340,1000,50],'BackgroundColor',[1 1 1],'FontSize',15,'FontName',...
    'Eras Demi ITC');
handles.title6=uicontrol(t1,'Style','edit','string',A(5)...
    ,'Position',[150,290,1000,50],'BackgroundColor',[1 1 1],'FontSize',15,'FontName',...
    'Eras Demi ITC');
handles.title7=uicontrol(t1,'Style','edit','string',A(6)...
    ,'Position',[150,240,1000,50],'BackgroundColor',[1 1 1],'FontSize',15,'FontName',...
    'Eras Demi ITC');
handles.title8=uicontrol(t1,'Style','edit','string',A(7)...
    ,'Position',[150,190,1000,50],'BackgroundColor',[1 1 1],'FontSize',15,'FontName',...
    'Eras Demi ITC');
handles.title9=uicontrol(t1,'Style','edit','string',A(8)...
    ,'Position',[150,140,1000,50],'BackgroundColor',[1 1 1],'FontSize',15,'FontName',...
    'Eras Demi ITC');
handles.title10=uicontrol(t1,'Style','edit','string',A(9)...
    ,'Position',[150,90,1000,50],'BackgroundColor',[1 1 1],'FontSize',15,'FontName',...
    'Eras Demi ITC');
%% Buttons for tab number 2
handles.Cake_push=uicontrol(t2,'Style','pushbutton','string','My Cake',...
    'Position',[15,500,120,50],'BackgroundColor',[1 1 1],'FontSize',12,'FontName',...
    'Eras Demi ITC','Callback',{@Cake_push_Callback});
handles.Candles_push=uicontrol(t2,'Style','pushbutton','string','Candles',...
    'Position',[15,440,120,50],'BackgroundColor',[1 1 1],'FontSize',12,'FontName',...
    'Eras Demi ITC','enable','off','Callback',{@Candles_push_Callback});
handles.Blow_off_push=uicontrol(t2,'Style','pushbutton','string','Blow off',...
    'Position',[15,380,120,50],'BackgroundColor',[1 1 1],'FontSize',12,'FontName',...
    'Eras Demi ITC','enable','off','Callback',{@Blow_off1_push_Callback});
%% Callbacks
function Cake_push_Callback(source,eventdata)
axes1=axes('Parent',t2);
[x_y y_y z_y]=cylinder(1,100);
z_y=z_y/3+0.5;
[x_r y_r z_r]=cylinder(1,100);
z_r=z_r/2;
surf(x_y,y_y,z_y,'facecolor',[1 1 0],'linestyle','none');  hold on
surf(x_r,y_r,z_r,'facecolor',[0.5 1 1],'linestyle','none'); grid off; axis on
fill(cos(linspace(0,2*pi,1000)),sin(linspace(0,2*pi,1000)),'y');
z_y_max=max(z_y);
z_y_max=z_y_max(1);
x_top=cos(linspace(0,2*pi,1000));
y_top=sin(linspace(0,2*pi,1000));
z_top=ones(size(x_top))*z_y_max;
plot3(x_top,y_top,z_top,'y');
fill3(x_top,y_top,z_top,'y');
x_taz=linspace(-0.65,0.65,1000);
y_taz=sin(20*x_taz)/10+0.6;
plot3(x_taz,y_taz,z_top,'r','linewidth',5)
x_taz=linspace(-0.95,0.95,1000);
y_taz=sin(20*x_taz)/10;
plot3(x_taz,y_taz,z_top,'r','linewidth',5)
x_taz=linspace(-0.65,0.65,1000);
y_taz=sin(20*x_taz)/10+0.6;
plot3(x_taz,-y_taz,z_top,'r','linewidth',5)
plot3(0.95*x_top,0.95*y_top,z_top,'r','linewidth',5)
% CHEERIES
[x_che y_che z_che]=sphere(100);
x_che=x_che*0.05;
y_che=y_che*0.05;
z_che=z_che*0.05;
for i=-0.8643:0.3138:0.7067
    surf(x_che+i,y_che-0.08,z_che+z_y_max,'facecolor',[0.6 0 0],'linestyle','none')
end
for i=-0.5385:0.3110:0.4
    surf(x_che+i,y_che-0.55,z_che+z_y_max,'facecolor',[0.6 0 0],'linestyle','none')
end
% Writing NOOSHIN
plot3([-0.6 -0.6]-0.1,[-0.4 -0.2],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([-0.4 -0.4]-0.1,[-0.4 -0.2],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([-0.6 -0.4]-0.1,[-0.2 -0.4],[z_y_max z_y_max],'m', 'linewidth',5)
plot3(0.1*cos((0:0.01:2*pi))-0.35,0.1*sin((0:0.01:2*pi))-0.3,...
    z_y_max*ones(629),'m', 'linewidth',5)
plot3(0.1*cos((0:0.01:2*pi))-0.12,0.1*sin((0:0.01:2*pi))-0.3,...
    z_y_max*ones(629),'m', 'linewidth',5)
plot3(0.1*cos((pi/4:0.01:pi))+0.12,0.1*sin((pi/4:0.01:pi))-0.3,...
    z_y_max*ones(236),'m', 'linewidth',5)
plot3(0.1*cos((pi+pi/4:0.01:2*pi))+0.12,-0.1*sin((pi/4:0.01:pi))-0.3,...
    z_y_max*ones(236),'m', 'linewidth',5)
plot3([0.37 0.37]-0.1,[-0.4 -0.2],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([0.47 0.47]-0.1,[-0.4 -0.2],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([0.37 0.47]-0.1,[-0.3 -0.3],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([0.55 0.55]-0.1,[-0.4 -0.2],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([0.63 0.63]-0.1,[-0.4 -0.2],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([0.83 0.83]-0.1,[-0.4 -0.2],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([0.63 0.83]-0.1,[-0.2 -0.4],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([0.12 0.32]-0.1,[-0.3 -0.3],[z_y_max z_y_max],'m', 'linewidth',5)
% Writing HAPPY
plot3([-0.6 -0.6]-0.1,[0.3 0.5],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([-0.5 -0.5]-0.1,[0.3 0.5],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([-0.6 -0.5]-0.1,[0.4 0.4],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([-0.45 -0.4]-0.1,[0.3 0.5],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([-0.4 -0.35]-0.1,[0.5 0.3],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([-0.42 -0.37]-0.1,[0.4 0.4],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([-0.28 -0.28]-0.1,[0.3 0.5],[z_y_max z_y_max],'m', 'linewidth',5)
plot3(0.05*cos(0:0.01:2*pi)-0.31,0.05*sin(0:0.01:2*pi)+0.45,...
    z_y_max*ones(629),'m', 'linewidth',5)
plot3([-0.12 -0.12]-0.1,[0.3 0.5],[z_y_max z_y_max],'m', 'linewidth',5)
plot3(0.05*cos(0:0.01:2*pi)-0.15,0.05*sin(0:0.01:2*pi)+0.45,...
    z_y_max*ones(629),'m', 'linewidth',5)
plot3([0.05 0.1]-0.1,[0.5 0.4],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([0.1 0.15]-0.1,[0.4 0.5],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([0.1 0.1]-0.1,[0.4 0.3],[z_y_max z_y_max],'m', 'linewidth',5)
% Writing BIRTHDAY
plot3([-0.1 -0.1]-0.2,[0.1 0.3],[z_y_max z_y_max],'m', 'linewidth',5)
plot3(0.05*cos(0:0.01:2*pi)-0.25,0.05*sin(0:0.01:2*pi)+0.275,...
    z_y_max*ones(629),'m', 'linewidth',5)
plot3(0.06*cos(0:0.01:2*pi)-0.24,0.06*sin(0:0.01:2*pi)+0.165,...
    z_y_max*ones(629),'m', 'linewidth',5)
plot3([0.07 0.07]-0.2,[0.1 0.3],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([0.15 0.15]-0.2,[0.1 0.3],[z_y_max z_y_max],'m', 'linewidth',5)
plot3(0.05*cos(0:0.01:2*pi),0.05*sin(0:0.01:2*pi)+0.275,...
    z_y_max*ones(629),'m', 'linewidth',5)
plot3([0.3 0.2]-0.2,[0.1 0.23],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([0.4 0.4]-0.2,[0.1 0.3],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([0.31 0.49]-0.2,[0.3 0.3],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([0.5 0.5]-0.2,[0.1 0.3],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([0.6 0.6]-0.2,[0.1 0.3],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([0.6 0.5]-0.2,[0.2 0.2],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([0.67 0.67]-0.2,[0.09 0.31],[z_y_max z_y_max],'m', 'linewidth',5)
plot3(0.05*cos(0:0.01:2*pi)+0.5,0.1*sin(0:0.01:2*pi)+0.21,...
    z_y_max*ones(629),'m', 'linewidth',5)
plot3([0.8 0.85]-0.2,[0.1 0.3],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([0.85 0.9]-0.2,[0.3 0.1],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([0.82 0.87]-0.2,[0.2 0.2],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([0.95 1]-0.2,[0.3 0.2],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([1 1.05]-0.2,[0.2 0.3],[z_y_max z_y_max],'m', 'linewidth',5)
plot3([1 1]-0.2,[0.1 0.2],[z_y_max z_y_max],'m', 'linewidth',5)
% CANDLES & FLAMES
for i=1:20
   [x_can y_can z_can]=cylinder(0.015,10);
   z_can=z_can/4+z_y_max;
   x_can=x_can+0.95*cos(i*pi/10);
   y_can=y_can+0.95*sin(i*pi/10);
   surf(x_can,y_can,z_can,'facecolor',[unifrnd(0.8,1,1) unifrnd(0,0.8,1) ...
      unifrnd(0,0.8,1)]); axis equal
end
set(handles.Candles_push,'enable','on');
set(handles.Blow_off_push,'enable','on');
set(handles.Cake_push,'enable','off');
msgbox('Press OK, then click the "Candles" button to turn them on.');
end
function Candles_push_Callback(source,eventdata)
    for i=1:20
      [x_f y_f z_f]=ellipsoid(0,0,0,0.01,0.01,0.08);
      x_f=x_f+0.95*cos(i*pi/10);
      y_f=y_f+0.95*sin(i*pi/10);
      AA=1.0833;
      z_f=z_f+AA;
      surf(x_f,y_f,z_f,'facecolor',[1 0.4 0.2],'linestyle','none'); axis equal
    end
    msgbox('Now the candles are on.')
end
function Blow_off1_push_Callback(source,eventdata)
        answer=questdlg('You must wish first, did you?','Question','Yes','No','No');
        if answer=='Yes'
            for i=1:20
      [x_f y_f z_f]=ellipsoid(0,0,0,0.01,0.01,0.08);
      x_f=x_f+0.95*cos(i*pi/10);
      y_f=y_f+0.95*sin(i*pi/10);
      AA=1.0833;
      z_f=z_f+AA;
      surf(x_f,y_f,z_f,'facecolor',[1 1 1],'linestyle','none'); axis equal
            end
            msgbox('The candles went off. HAPPY BIRTHDAY again ! ! !.');
        end
end
%% Buttons for tab number 3
A={'Good luck & Bye !!!)'...
 'Esmaeel Kakavand',...
'Sunday - - - 1391/Farvardin/6 - - - 2012/Mars/26'};
handles.title12=uicontrol(t3,'Style','edit','string',A(1)...
    ,'Position',[150,490,1000,50],'BackgroundColor',[1 1 1],'FontSize',15,'FontName',...
    'Eras Demi ITC');
handles.title13=uicontrol(t3,'Style','edit','string',A(2)...
    ,'Position',[150,440,1000,50],'BackgroundColor',[1 1 1],'FontSize',15,'FontName',...
    'Eras Demi ITC');
handles.title14=uicontrol(t3,'Style','edit','string',A(3)...
    ,'Position',[150,390,1000,50],'BackgroundColor',[1 1 1],'FontSize',15,'FontName',...
    'Eras Demi ITC');
handles.quit_push=uicontrol(t3,'Style','pushbutton','string','Quit'...
    ,'Position',[600,320,100,50],'BackgroundColor',[1 1 1],'FontSize',15,'FontName',...
    'Eras Demi ITC','BackgroundColor','y','Callback',{@Quit_push_Callback});
    function Quit_push_Callback(source,eventdata)
        close;
    end
end