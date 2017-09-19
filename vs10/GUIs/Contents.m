% GUIs - building and using Graphical User Interfaces
%   
%  Creating a GUI
%     newGUI       - open GUI figure
%     GUIsettings  - settings for GUIs
%
%  GUI components
%     GUIpanel     - define, but do not draw, uipanel containing param queries
%     ParamQuery   - constructor for ParamQuery objects
%     ActionButton - constructor for ActionButton objects
%     toggle       - constructor for toggle object
%     messenger    - constructor for Messenger objects
%     AxesDisplay  - constructor for AxesDisplay objects
%     cycleList    - constructor for CycleList objects
%     pulldownmenu - constructor for PulldownMenu objects
%     GUIpiece     - construct GUIpiece object
%
%  Reading and writing data
%     GUIval       - evaluate user input from GUI
%     GUImessage   - display general message in GUI
%     GUIgrab      - get strings from GUI edits and toggles without interpreting them
%     GUIfill      - restore stored string values of paramqueries of a GUI
%     getGUIdata   - get struct field of userdata property of a graphics object(from getguidata)
%     rmGUIdata    - remove struct field of userdata property of a graphics object(from rmguidata)
%     setGUIdata   - set struct field of userdata property of a graphics object(from setguidata)
%
%  Retrieving GUI components
%     gcg          - handle to current GUI
%     findGUIobj   - find named messenger object in GUI
%     GUIaxes      - handle of axes in GUI
%     GUImessenger - find named messenger object in GUI
%     isCycleList  - true for cycle list
%
%  See also GUIutils.


