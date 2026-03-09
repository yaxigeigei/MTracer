# MSessionExplorer (MATLAB Version)

This is the MATLAB implementation of MSessionExplorer. Note that the MATLAB version is completely different from the Python version.

## Installation

1. Clone or download this repository
2. Add the `MATLAB` subfolder to your MATLAB path:
```matlab
addpath('path/to/MSessionExplorer/MATLAB');
```

## Overview

MSessionExplorer is designed to handle three types of data tables, where each row represents an epoch (trial/session):

1. **Time Series Tables**
   - Store continuous time series data
   - First column must be named 'time' and contains cell arrays of timestamps
   - Additional columns contain cell arrays of time series values
   - Each cell contains the data for one epoch
   - Example:
   ```matlab
   time = {[0:0.1:1]'; [0:0.1:2]'};
   signal = {sin(0:0.1:1)'; cos(0:0.1:2)'};
   timeseries_table = cell2table([time, signal], 'VariableNames', {'time', 'signal'});
   ```

2. **Event Times Tables**
   - Store timestamps of discrete events
   - Each column represents a different event type
   - Values are stored in cell arrays, where each cell contains the timestamps for one epoch
   - Values can be simple timestamps or `MSessionExplorer.Event` objects that carry additional metadata
   - Example: 
   ```matlab
   event_A = {[1.0, 2.0, 3.0]'; [1.5, 2.5]'};
   event_B = {[1.2, 2.2]'; NaN};
   event_times_table = cell2table([event_A, event_B], 'VariableNames', {'event_A', 'event_B'});
   ```

3. **Event Values Tables**
   - Store event-related data that doesn't have a time component
   - Each row corresponds to an epoch
   - Values can be any MATLAB data type
   - Example:
   ```matlab
   trial_num = {1; 2};
   condition = {'A'; 'B'};
   success = {true; false};
   event_values_table = cell2table([trial_num, condition, success], ...
       'VariableNames', {'trial_num', 'condition', 'success'});
   ```

## Key Features

1. **Data Management**
   - Store and organize multiple types of data tables
   - Track reference times for temporal alignment
   - Store arbitrary metadata in the `userData` property
   - Support for epochs/trials with automatic synchronization

2. **Data Operations**
   - Time alignment and morphing
   - Data slicing and resampling
   - Epoch/trial management (sorting, splitting, merging)
   - Event time manipulation with the `MSessionExplorer.Event` class

3. **Event Class**
   - Special class for handling event timestamps with metadata
   - Supports time operations (addition, subtraction, scaling)
   - Can store additional timestamps and arbitrary data
   - Useful for complex event-based analysis

## Documentation and Examples

You can open the documentation pages for the `MSessionExplorer` class and the `MSessionExplorer.Event` class by executing the following lines in the MATLAB command window:
```matlab
doc MSessionExplorer
doc MSessionExplorer.Event
```

Extensive examples can be accessed by executing the following lines in the MATLAB command window:
```matlab
MSessionExplorer.Examples.MakeEventTimesTable
MSessionExplorer.Examples.MakeTimeSeriesTable
MSessionExplorer.Examples.Duplicate
MSessionExplorer.Examples.Merge
MSessionExplorer.Examples.Peek
MSessionExplorer.Examples.SetTable
MSessionExplorer.Examples.SetReferenceTime
MSessionExplorer.Examples.SetColumn
MSessionExplorer.Examples.AlignTime
MSessionExplorer.Examples.SliceSession
MSessionExplorer.Examples.SliceEventTimes
MSessionExplorer.Examples.SliceTimeSeries
MSessionExplorer.Examples.ResampleEventTimes
MSessionExplorer.Examples.ResampleTimeSeries
```