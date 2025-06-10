classdef MLing
    % MLing is...
    
    methods(Static)
        function [k1, k2, alignInfo] = FindAlignedWords(s1, s2)
            % Globally align two sentences at the level of letters, then find pairs of corresponding words.
            % 
            %   [k1, k2, alignInfo] = MLing.FindAlignedWords(s1, s2)
            % 
            % Inputs
            %   s1, s2          Two sequences of words (e.g. sentences). If an input is a vector of 
            %                   strings, elements will be joined together by spaces into a sentence.
            % Outputs
            %   k1, k2          Indices of matching words in s1 and s2, respectively. Namely, the 
            %                   k1(i)-th word in s1 matches the k2(i)-th word in s2.
            %   alignInfo       A struct storing the alignment result in the following fields.
            %                   s1, s2                  The two input sequences. Concatenated if any 
            %                                           of the input is a vector of words.
            %                   aa, steps, tscore       Outputs from MLing.SeqAlign function. See 
            %                                           its documentation for details.
            % 
            
            ss = {s1, s2};
            wdInd = cell(size(ss));
            for i = 1 : numel(ss)
                s = ss{i};
                
                % Find which word (by index) each character belongs to
                s = strtrim(s);
                if ~ischar(s) && numel(s) > 1
                    w = s;
                    s = strjoin(w, ' ');
                else
                    w = strsplit(s, ' ');
                end
                wLen = strlength(w) + 1; % plus one to include a trailing space
                wInd = repelem(1:numel(w), wLen);
                wInd(end) = []; % remove index of the trailing space of the last word
                
                ss{i} = s;
                wdInd{i} = wInd;
            end
            
            % Global alignment
            [aa, steps, tscore] = MLing.SeqAlign(char(ss{1}), char(ss{2}));
            
            % Find each character's word index in the aligned strings
            chInd = cumsum(steps, 2);
            chInd(~chInd) = 1; % avoid zero indexing
            wdInd = [wdInd{1}(chInd(1,:)); wdInd{2}(chInd(2,:))];
            
            % Check word pairs
            isGap = any(~steps, 1);
            wdIndNoGap = wdInd(:, ~isGap);
            prInd = unique(wdIndNoGap', 'rows', 'stable');
            
            % Output
            k1 = prInd(:,1);
            k2 = prInd(:,2);
            alignInfo.s1 = ss{1};
            alignInfo.s2 = ss{2};
            alignInfo.aa = aa;
            alignInfo.steps = steps;
            alignInfo.tscore = tscore;
        end
        
        function [k1, k2, alignInfo] = FindAlignedTokens(s1, s2)
            % Globally align two sequences of tokens, and find pairs of corresponding ones.
            % 
            %   [k1, k2, alignInfo] = MLing.FindAlignedTokens(s1, s2)
            % 
            % Inputs
            %   s1, s2          Two sequences of tokens (e.g. phonemes), each is a vector of string
            %                   objects or a cell vector of char strings.
            % Outputs
            %   k1, k2          Indices of matching tokens in s1 and s2, respectively. Namely, the 
            %                   k1(i)-th token in s1 matches the k2(i)-th token in s2.
            %   alignInfo       A struct storing the alignment result in the following fields.
            %                   s1, s2                  The two input sequences.
            %                   aa, steps, tscore       Outputs from MLing.SeqAlign function. See 
            %                                           its documentation for details.
            % 
            % See also MLing.SeqAlign
            
            % Global alignment
            s1 = string(s1);
            s2 = string(s2);
            [aa, steps, tscore] = MLing.SeqAlign(s1, s2);
            
            % Get token indices
            tkInd = cumsum(steps, 2);
            
            % Find token pairs
            isGap = any(~steps, 1);
            tkInd = tkInd(:, ~isGap);
            
            % Output
            k1 = tkInd(1,:)';
            k2 = tkInd(2,:)';
            alignInfo.s1 = s1;
            alignInfo.s2 = s2;
            alignInfo.aa = aa;
            alignInfo.steps = steps;
            alignInfo.tscore = tscore;
        end
        
        function id = FindTimitID(srcDir, varargin)
            % Find TIMIT stim IDs
            % 
            %   id = FindTimitID(srcDir)
            %   id = FindTimitID(srcDir, idPattern)
            % 
            % Inputs
            %   srcDir      The source directory path of TIMIT feature/label files.
            %   idPattern   A string of text pattern to match for. Default is '*_*.wav'.
            % Output
            %   id          A cell or string array of TIMIT stim IDs.
            % 
            
            assert(exist(srcDir, 'dir'), "The source directory does not exist.");
            
            p = inputParser();
            p.addOptional('pattern', '*_*.wav', @(x) isstring(x) || ischar(x) || iscellstr(x));
            p.parse(varargin{:});
            pattern = p.Results.pattern;
            
            % Search for all sentences based on wav file
            search = MBrowse.Dir2Table(fullfile(srcDir, pattern));
            [~, id] = fileparts(search.name);
            
            id = cellstr(id);
            
            % Exclude bad sentences
            badIds = {'mdwh0_si1925'};
            isBad = ismember(id, badIds);
            if any(isBad)
                fprintf("Exclude %s due to incomplete labels\n", id{isBad});
                id(isBad) = [];
            end
        end
        
        function [wf, t] = ReadTimitWaveform(srcDir, id)
            % Read TIMIT audio waveform
            %
            %   [wf, t] = MLing.ReadTimitWaveform(timitDir, id)
            % 
            % Inputs
            %   srcDir      The source directory path of TIMIT feature/label files.
            %   id          A cell or string array of TIMIT stim IDs.
            % Outputs
            %   wf          A cell array of audio waveform.
            %   t           Timestamps associated with wf in seconds.
            % 
            assert(exist(srcDir, 'dir'), "The source directory does not exist.");
            id = cellstr(id);
            wf = cell(size(id));
            t = wf;
            fs = 16e3; % hardcode the TIMIT sampling frequency
            for i = 1 : numel(wf)
                p = fullfile(srcDir, id{i} + ".wav");
                wf{i} = audioread(p);
                t{i} = (1:length(wf{i}))' / fs;
            end
        end
        
        function tg = ReadTimitFeatures(srcDir, varargin)
            % Construct textgrid compatible struct from wrd, phn
            %
            %   tg = MLing.ReadTimitFeatures(timitDir, id)
            % 
            % Inputs
            %   srcDir      The source directory path of TIMIT feature/label files.
            %   id          A list of TIMIT stim IDs.
            % Output
            %   tg          An array of structs compatible to use with the TextGrid MATLAB package.
            %               Currently it only contains words and phones level.
            % 
            % TODO          Add syllable tier to tg.
            %               Read quantity features such as pitch and formants.
            % 
            
            p = inputParser();
            p.addOptional('id', [], @(x) isstring(x) || ischar(x) || iscellstr(x));
            p.addParameter('UniformOutput', true, @islogical);
            p.parse(varargin{:});
            id = p.Results.id;
            isUni = p.Results.UniformOutput;
            
            id = cellstr(id);
            tg = cell(size(id));
            
            for i = 1 : numel(id)
                % Read wrd and phn files as tables
                wrdFile = fullfile(srcDir, id{i} + ".wrd");
                phnFile = fullfile(srcDir, id{i} + ".phn");
                if ~exist(wrdFile, 'file') || ~exist(phnFile, 'file')
                    warning('Cannot find wrd and/or phn files with stim ID: ''%s''', id{i});
                    continue
                end
                wrdTb = readtable(wrdFile, 'FileType', 'text');
                phnTb = readtable(phnFile, 'FileType', 'text');
                
                % Remove silent phones
                phnTb(strcmp(phnTb.(3), 'h#'), :) = [];
                phnTb.(3) = upper(phnTb.(3));
                
                % Hardcode the audio sampling rate of TIMIT
                fs = 16000;
                
                % Construct TextGrid compatible struct
                w = struct;
                w.name = 'words';
                w.type = 'interval';
                w.T1 = wrdTb.(1) / fs;
                w.T2 = wrdTb.(2) / fs;
                w.Label = wrdTb.(3)';
                
                p = struct;
                p.name = 'phones';
                p.type = 'interval';
                p.T1 = phnTb.(1) / fs;
                p.T2 = phnTb.(2) / fs;
                p.Label = phnTb.(3)';
                
                s = struct;
                s.tier = {w, p};
                s.tmin = w.T1(1);
                s.tmax = w.T2(end);
                
                tg{i} = s;
            end
            
            % Convert to struct array or remain as cell array
            if isUni
                tg = cat(1, tg{:});
            end
        end
        
        function labels = ARPA2IPA(labels)
            % Map ARPAbet symbols to IPA symbols
            % 
            %   labels = MLing.ARPA2IPA(labels)
            % 
            
            dict = MLing.dict(:, [1 2]);
            
            dtype = class(labels);
            
            labels = cellstr(labels);
            for k = 1 : numel(labels)
                L = upper(labels{k});
                L(L>='0' & L<='9') = []; % remove the trailing digit in phonemes
                I = find(strcmp(L, dict(:,1)), 1);
                if ~isempty(I)
                    labels{k} = dict{I,2};
                end
            end
            
            if dtype == "char"
                labels = labels{1};
            elseif dtype == "string"
                labels = string(labels);
            end
        end
        
    end
    
    properties(Constant)
        % Symbol conversion lookup
        dict = { % columns: ARPAbet, IPA
            % Standard
            'AA',   'ɑ';
            'AE',   'æ';
            'AH',   'ʌ';
            'AO',   'ɔ';
            'AW',   'aʊ';
            'AX',   'ə';
            'AXR',  'ɚ';
            'AY',   'aɪ';
            'B',    'b';
            'CH',   'tʃ';
            'D',    'd';
            'DH',   'ð';
            'DX',   'ɾ';
            'EL',   'l̩';
            'EM',   'm̩';
            'EN',   'n̩';
            'EH',   'ɛ';
            'ER',   'ɝ';
            'EY',   'eɪ';
            'F',    'f';
            'G',    'ɡ';
            'HH',   'h';
            'IH',   'ɪ';
            'IX',   'ɨ';
            'IY',   'i';
            'JH',   'dʒ';
            'K',    'k';
            'L',    'l';
            'M',    'm';
            'N',    'n';
            'NG',   'ŋ';
            'NX',   'ɾ̃';
            'OW',   'oʊ';
            'OY',   'ɔɪ';
            'P',    'p';
            'Q',    'ʔ';
            'R',    'ɹ';
            'S',    's';
            'SH',   'ʃ';
            'T',    't';
            'TH',   'θ';
            'UH',   'ʊ';
            'UW',   'u';
            'UX',   'ʉ';
            'V',    'v';
            'W',    'w';
            'Y',    'j';
            'Z',    'z';
            'ZH',   'ʒ';
            % TIMIT specific
            'AX-H', 'ə̥';
            'BCL',  'b̚';
            'DCL',  'd̚';
            'ENG',  'ŋ̍';
            'GCL',  'ɡ̚';
            'HV',   'ɦ';
            'KCL',  'k̚';
            'PCL',  'p̚';
            'TCL',  't̚';
            'PAU',  '';
            'EPI',  '';
            'H#',   '';
            };
    end
    
end

