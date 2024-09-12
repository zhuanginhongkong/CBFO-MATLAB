% This Class aims to produce GUI for CBFO
classdef GUITopo
  methods(Static)
    % Create Handle for all Interface Window
    function  Handle = GUIDemo(obj)
      % Create Window and
      Handle.Demoobj = obj ;
      Handle.Window = figure(3) ;
      Handle.Window.Color = [1 1 1]*0.96 ; Handle.Window.MenuBar = 'none' ;
      Handle.Window.Name = 'Result Display' ;  Handle.Window.NumberTitle = 'off' ;
      Handle.Window.Units = 'normalized';
      Handle.Window.Position =[0.63 0.05  0.35 0.9] ;
      Handle.TEXT{1} = uicontrol(Handle.Window,'style','text','String','First Remeshing',...
        'units','normalized','fontsize',16) ;
      Handle.TEXT{1}.Position = [0.05 0.95 0.9 0.03] ;
      Handle.TEXT{1}.BackgroundColor = [1 1 1]*0.96 ;
      Handle.TEXT{1}.HorizontalAlignment = 'left' ;
      Handle.TEXT{2} = uicontrol(Handle.Window,'style','text','String','Second Remeshing',...
        'units','normalized','fontsize',16) ;
      Handle.TEXT{2} .Position = [0.05 0.65 0.9 0.03] ; Handle.TEXT{2} .BackgroundColor = [1 1 1]*0.96 ;
      Handle.TEXT{2} .HorizontalAlignment = 'left' ;
      % Create Subfigure
      Handle.Subfigure{1} = axes(Handle.Window,'units','normalized','position',[0.09 0.71  0.75 0.22 ]) ;
      axis(Handle.Subfigure{1},'off')
      Handle.Subfigure{2} = axes(Handle.Window,'units','normalized','position',[0.09 0.41  0.75 0.22 ]) ;
      axis(Handle.Subfigure{2},'off')
      % Create Slider
      Handle.Slider{1} = uicontrol(Handle.Window,'style','slider','units','normalized',...
        'position',[0.88 0.71 0.04 0.23]) ;
      Handle.Slider{2} = uicontrol(Handle.Window,'style','slider','units','normalized',...
        'position',[0.88 0.41 0.04 0.23]) ;
      Handle.Slider{1}.Callback = @(HObject, Eventdata) ...
        GUITopo.ResultDemo(HObject, Eventdata, Handle) ;
      Handle.Slider{2}.Callback = @(HObject, Eventdata) ...
        GUITopo.ResultDemo(HObject, Eventdata, Handle) ;
    end
    %% This is the Callback function to demo the topology result
    function HObject = ResultDemo(HObject, Eventdata,Handle)
      % Link the static function
      Display = @ GUITopo.IterationDemo ;
      % Get the Iteration times for Initial meshing and First body-fitted Remeshing 
      Num1 = size(Handle.Demoobj.Output.FirstMeshing.Xphy, 2)  ;
      Num2 = size(Handle.Demoobj.Output.SecondMeshing.Xphy, 2) ;
      FitstIteration = max(round(Handle.Slider{1}.Value*Num1), 1) ;
      SecondIteration = max(round(Handle.Slider{2}.Value *Num2), 1) ;
      % Change the static text
      Handle.TEXT{1}.String = [num2str(FitstIteration), ' Iterations after initial meshing'  ] ;
      Handle.TEXT{2}.String = [num2str(SecondIteration), ' Iterations after first body-fitted remeshing'   ] ;
      % Demo the Iteration Process
      pFirst = Handle.Demoobj.Output.FirstMeshing.p ;
      tFirst = Handle.Demoobj.Output.FirstMeshing.t ;
      XFirst = Handle.Demoobj.Output.FirstMeshing.Xphy{FitstIteration} ;
      pSecond = Handle.Demoobj.Output.SecondMeshing.p ;
      tSecond = Handle.Demoobj.Output.SecondMeshing.t ;
      XSecond = Handle.Demoobj.Output.SecondMeshing.Xphy{SecondIteration} ;
      Display(pFirst,tFirst,XFirst, Handle.Subfigure{1}) ;
      Display(pSecond,tSecond,XSecond, Handle.Subfigure{2}) ;
    end
    % This function is to demo the optimization progress before the first and second remeshing
    function IterationDemo(p, t, xnew, Destination)
      colormap(Destination,'summer') ;
      hold(Destination,'off')
      patch( Destination, 'Faces',t, 'Vertices', p, 'FaceVertexCData', ...
        xnew, 'FaceColor', 'flat') ;
      axis(Destination,'equal') ; axis(Destination,'off')
    end
  end
end