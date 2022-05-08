classdef OperationManager < handle
    %OperationManager Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        stateTb table;
        actionTb table;
        isVerbose = true;
    end
    
    methods
        function this = OperationManager(stateData)
            % OperationManager Construct an instance of this class
            
            % Initilize state table
            tb = table;
            tb.idx = 1;             % state index
            tb.data = stateData;    % a struct of state data
            this.stateTb = tb;
            
            % Initialize action table
            tb = table;
            tb.idx = 1;             % action index
            tb.head = 1;            % head location in action index
            tb.fromState = 1;       % action input state index
            tb.toState = 1;         % action output state index
            tb.type = "init";     % type of the action: 'init', 'op', 'undo', or 'redo'
            this.actionTb = tb;
            
            if this.isVerbose
                disp(this.actionTb);
            end
        end
        
        function AddOperation(this, stateData)
            % Add an operation with new stateData
            
            % Append new state to the state table
            tb = this.stateTb;
            idx = height(tb) + 1;
            this.stateTb = [tb; {idx, stateData}];
            
            % Append action to the action table
            tb = this.actionTb;
            idx = height(tb) + 1;
            head = idx;
            fromState = tb.toState(idx-1);
            toState = height(this.stateTb);
            type = "op";
            this.actionTb = [tb; {idx, head, fromState, toState, type}];
            
            if this.isVerbose
                disp(stateData);
                disp(this.actionTb);
            end
        end
        
        function [fromStateData, toStateData] = Undo(this)
            % Add an undo action and return the input output state data
            
            tb = this.actionTb;
            
            % Check if the head is already the earliest
            if tb.head(end) == 1
                fromStateData = [];
                toStateData = [];
                if this.isVerbose
                    disp('No further action can be undone.');
                end
                return
            end
            
            % Append action to the action table
            idx = height(tb) + 1;
            head = tb.head(end) - 1;                % decrement the head from last action
            fromState = tb.toState(tb.head(end));   % reverse states of the action at last head
            toState = tb.fromState(tb.head(end));   % reverse states of the action at last head
            type = "undo";
            this.actionTb = [tb; {idx, head, fromState, toState, type}];
            
            % Get state data
            fromStateData = this.stateTb.data(fromState);
            toStateData = this.stateTb.data(toState);
            
            if this.isVerbose
                disp(this.actionTb);
            end
        end
        
        function [fromStateData, toStateData] = Redo(this)
            % Add an redo action and return the input output state data
            
            tb = this.actionTb;
            
            % Check if the head reached the last op
            if tb.head(end) == max(tb.head(tb.type=="op"))
                fromStateData = [];
                toStateData = [];
                if this.isVerbose
                    disp('No further action can be redone.');
                end
                return
            end
            
            % Append action to the action table
            idx = height(tb) + 1;
            head = tb.head(end) + 1;        % increment the head from last action
            fromState = tb.fromState(head); % repeat states of the action at head
            toState = tb.toState(head);     % repeat states of the action at head
            type = "redo";
            this.actionTb = [tb; {idx, head, fromState, toState, type}];
            
            % Get state data
            fromStateData = this.stateTb.data(fromState);
            toStateData = this.stateTb.data(toState);
            
            if this.isVerbose
                disp(this.actionTb);
            end
        end
    end
end

