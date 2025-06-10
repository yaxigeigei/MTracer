# MSessionExplorer
A tool that helps munging heterogenous timeseries and event-based data for your analysis.

## Versions
**This README introduces the Python version of MSessionExplorer.**  
For more detailed information about the MATLAB version, please see the README in the `@MATLAB` folder.

## Installation

### Python Version
```bash
pip install MSessionExplorer
```

## Data Structures

### Data Object Properties
The `Data` class provides several useful properties:
- `table_names`: List of all table names in the session
- `table_types`: List of table types ('timeseries', 'interval', or 'timestamp')
- `dod` (DataFrame of DataFrames): A summary DataFrame showing table information including:
  - Table name
  - Table type
  - Number of rows and columns
  - Column names

### Table Types and Required Columns

1. **Timeseries Tables**
   - Must have a 'time' column as the first column
   - Additional columns can contain any numeric data
   - Example structure:
   ```python
   pd.DataFrame({
       'time': [0.0, 0.1, 0.2, ...],  # Required first column
       'series1': [1.0, 2.0, 3.0, ...],  # Optional data columns
       'series2': [4.0, 5.0, 6.0, ...]
   })
   ```

2. **Interval Tables**
   - Must have 'start_time' and 'stop_time' as the first two columns
   - Additional columns can contain any data type
   - Example structure:
   ```python
   pd.DataFrame({
       'start_time': [1.0, 3.0, 5.0],  # Required first column
       'stop_time': [2.0, 4.0, 6.0],   # Required second column
       'event_type': ['A', 'B', 'C'],  # Optional data columns
       'amplitude': [1.5, 2.0, 1.8]
   })
   ```

3. **Timestamp Tables**
   - Each column represents a different event type
   - Values are arrays of timestamps
   - Example structure:
   ```python
   pd.DataFrame({
       'event_A': [[1.0, 3.0, 5.0]],  # Each cell contains an array of timestamps
       'event_B': [[2.0, 4.0, 6.0]]
   })
   ```

### Metadata
The `metadata` property is a dictionary that can store any additional information about the session:
```python
se.metadata['session_id'] = 'session_001'
se.metadata['experimenter'] = 'John Doe'
se.metadata['notes'] = 'First recording session'
```

## Usage

### Data Class
The `Data` class helps manage session data in timeseries or event-based tables and metadata dictionaries.  
Here's a simple example:

```python
from MSessionExplorer import Data

# Create a new Data object
se = Data()

# Add a timeseries table
import pandas as pd
import numpy as np

# Create a sample timeseries
time = np.arange(0, 10, 0.1)
values = np.sin(time)
timeseries_df = pd.DataFrame({
    'time': time,
    'sine': values
})
se.set_table('sine_wave', timeseries_df)

# Add an interval table
intervals_df = pd.DataFrame({
    'start_time': [1.0, 3.0, 5.0],
    'stop_time': [2.0, 4.0, 6.0],
    'event_type': ['A', 'B', 'C']
})
se.set_table('events', intervals_df)

# View the data
se.peek()
```

### Slicer Class
The `Slicer` class helps extract data within specific time intervals. Here's an example:

```python
from MSessionExplorer import Slicer

# Create a slicer
slicer = Slicer(target_se=se)

# Create intervals that will be used to slice other data
slicer_intervals = pd.DataFrame({
    'start_time': [0.5, 2.5],
    'stop_time': [1.5, 3.5]
})
slicer.set_table('slicer', slicer_intervals)

# Slice the timeseries data
sliced_data = slicer.slice(
    slicer_table='slicer',
    target_table='sine_wave',
    sample_rate=10  # Optional: resample to 10 Hz
)

# Slice the interval data
sliced_events = slicer.slice(
    slicer_table='slicer',
    target_table='events',
    interval_include='either'  # Include events that overlap with the interval
)
```