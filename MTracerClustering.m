classdef MTracerClustering
    %MTRACERCLUSTERING Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        vm MTracerVM;
        kr NP.KilosortResult;
        
        clusUiTb table;                 % cluster table displayed in the GUI
        clusHisTb table;                % like kr.clusTb but has an extra spkInd column, and never delete old clusters
        spkHisTb table;                 % each column contains clusId of a step
    end
    
    methods
        function this = MTracerClustering(vm, kr)
            % 
            this.vm = vm;
            this.kr = kr;
        end
        
        function Merge(this, cid)
            % 
            
            
        end
        
        function Cut(this, cid)
            % 
            
            
        end
        
        function LogHistory(this)
            % 
            
            
        end
        
        function Undo(this)
            % 
            
            
        end
        
        function Redo(this)
            % 
            
            
        end
        
        function SetPoint(this)
            % 
            
            
        end
        
        function RemovePoint(this)
            % 
            
            
        end
        
        function PlotPolygon(this)
            % 
            
            
        end
    end
end

