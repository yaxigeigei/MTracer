import numpy as np
import pandas as pd
from pathlib import Path
import dill
import copy
from IPython.display import display


class Data:
    """
    A class for managing session data with tables and metadata.
    """
    
    @property
    def dod(self):
        return pd.DataFrame({
            'table_name': self.table_names,
            'table_type': self.table_types,
            'num_rows': [len(self.get_table(table_name)) for table_name in self.table_names],
            'num_cols': [len(self.get_table(table_name).columns) for table_name in self.table_names],
            'columns': [self.get_table(table_name).columns.tolist() for table_name in self.table_names]
        })
    
    @property
    def table_names(self):
        return list(self._dod.keys())
    
    @property
    def table_types(self):
        return [self.get_table_type(table_name) for table_name in self.table_names]

    def __init__(self):
        """Initialize an empty MSessionExplorer.Data object."""
        self._dod = {}
        self.metadata = {}
    
    def peek(self):
        """
        Peek at the first few rows of each table.
        """
        # Display summary of tables
        print('Tables:')
        display(self.dod)

        # Display metadata
        if not self.metadata:
            print('No metadata available')
            return
        
        for key, value in self.metadata.items():
            print(f"{key}:")
            display(value)

    def set_table(self, table_name, table_data):
        """
        Add or replace a dataframe in the session explorer.
        
        Parameters
        ----------
        table_name : str
            Name of the table to add or replace
        table_data : pandas.DataFrame
            DataFrame to store
        """
        self._dod[table_name] = table_data
    
    def get_table(self, *table_names):
        """
        Retrieve one or more dataframes from the session explorer.
        
        Parameters
        ----------
        *table_names : str
            Names of the tables to retrieve
            
        Returns
        -------
        pandas.DataFrame
            The requested dataframes
        """
        if not table_names:
            raise ValueError("No table names provided")
        elif len(table_names) == 1:
            return self._dod[table_names[0]].copy()
        else:
            return [self._dod[table_name].copy() for table_name in table_names]

    def get_table_type(self, *table_names):
        """
        Get the type of a table by name.
        """
        if not table_names:
            raise ValueError("No table names provided")
        
        table_types = []
        for table_name in table_names:
            df = self.get_table(table_name)
            table_types.append(self._get_table_type(df))
        
        if len(table_types) == 1:
            return table_types[0]
        else:
            return table_types
    
    def _get_table_type(self, df):
        """
        Get the type of a dataframe.
        """
        if df.columns[0] == 'time':
            return 'timeseries'
        elif df.columns[0] == 'start_time' and df.columns[1] == 'stop_time':
            return 'interval'
        else:
            return 'timestamp'

    def set_column(self, table_name, column_names, column_data):
        """
        Add or replace a column in a table.
        
        Parameters
        ----------
        table_name : str
            Name of the table to modify
        column_names : str, list of str, or int/list of int
            Name(s) of the column(s) to add, replace, or delete
        column_data : object or None
            Data to set for the column(s). If None, the column(s) will be deleted.
        """
        if not isinstance(column_names, (list, tuple)):
            column_names = [column_names]
        
        for idx, col in enumerate(column_names):
            if column_data is None:
                del self._dod[table_name][col]
            elif len(column_names) == 1:
                self._dod[table_name][col] = column_data
            else:
                self._dod[table_name][col] = column_data[idx]

    def get_column(self, table_name, column_names):
        """
        Get a column from a table.
        """
        if isinstance(column_names, (list, tuple)):
            return [self._dod[table_name][col] for col in column_names]
        else:
            return self._dod[table_name][column_names]

    def remove_table(self, *table_names):
        """
        Remove one or more tables from the session explorer.
        
        Parameters
        ----------
        *table_names : str
            Names of the tables to remove
            
        Returns
        -------
        list of pandas.DataFrame
            The removed dataframes in the same order as table_names
        """
        removed_tables = []
        for table_name in table_names:
            if table_name in self._dod:
                removed_tables.append(self._dod[table_name])
                del self._dod[table_name]
            else:
                print(f"Table '{table_name}' not found")
                removed_tables.append(None)
        return removed_tables

    def copy(self):
        """
        Create a deep copy of the MSessionExplorer.Data object.
        """
        return copy.deepcopy(self)

    def save(self, filepath):
        """
        Save the MSessionExplorer.Data object to a file.
        
        Parameters
        ----------
        filepath : str or Path
            Path to save the file
        """
        filepath = Path(filepath)
        with open(filepath, 'wb') as f:
            dill.dump(self, f)

    def morph_time(self, interpolant, decimals=6):
        """
        Morph the time-related columns of all tables.

        Parameters
        ----------
        interpolant : function
            A function that takes a time and returns a new time
        decimals : int, optional
            The number of decimal places to round to
        """

        def interpolate(x, n=decimals):
            return np.round(interpolant(x), n)

        for table_name in self.table_names:
            table = self._dod[table_name]
            table_type = self._get_table_type(table)
            if table_type == 'timeseries':
                self._dod[table_name]['time'] = interpolate(self._dod[table_name]['time'])
            elif table_type == 'interval':
                self._dod[table_name]['start_time'] = interpolate(self._dod[table_name]['start_time'])
                self._dod[table_name]['stop_time'] = interpolate(self._dod[table_name]['stop_time'])
            elif table_type == 'timestamp':
                for col in table.columns:
                    self._dod[table_name][col].iloc[0] = interpolate(self._dod[table_name][col].iloc[0])
            print(f"Morphed time for {table_type} table '{table_name}'")
    
    def resample_timeseries(self, table, new_time=None, sample_rate=None, inplace=False):
        """
        Resample a timeseries table.

        Parameters
        ----------
        table : str or pandas.DataFrame
            The table to resample
        new_time : array-like, optional
            The new timestamps
        sample_rate : float, optional
            The new sample rate
        inplace : bool, optional
            Whether to modify the table in place

        Returns
        -------
        pandas.DataFrame
            The resampled table
        """
        # Check table
        if isinstance(table, str):
            table_name = table
            table = self._dod[table]
        else:
            table_name = None
        assert isinstance(table, pd.DataFrame), "Table must be a pandas DataFrame"
        assert self._get_table_type(table) == 'timeseries', "Table must be of timeseries type"
        
        # Check resampling parameters
        if new_time is None and sample_rate is None:
            raise ValueError("Provide either new_time or sample_rate")
        elif new_time is not None and sample_rate is not None:
            raise ValueError("Provide either new_time or sample_rate, not both")
        
        # Get new timestamps
        old_time = table['time'].values
        if new_time is None:
            new_time = np.arange(old_time[0], old_time[-1], 1.0/sample_rate)
        
        # Resample timeseries
        new_table = {'time': new_time}
        for col in table.columns[1:]:
            old_data = table[col].values
            new_data = np.interp(new_time, old_time, old_data)
            new_table[col] = new_data
        new_table = pd.DataFrame(new_table)

        # Return new table
        if table_name is not None and inplace:
            self._dod[table_name] = new_table
        else:
            return new_table


class Slicer(Data):
    """
    A class for slicing MSessionExplorer.Data objects.
    """

    @property
    def target_se(self):
        return self._target_se
    
    @target_se.setter 
    def target_se(self, se):
        assert isinstance(se, Data) or se is None
        self._target_se = se

    def __init__(self, target_se=None):
        # Set the default target MSessionExplorer.Data object
        self.target_se = target_se
        
    def set_table(self, table_name, table_data):
        """
        Add or replace a dataframe in the session explorer.
        """
        if not self._get_table_type(table_data) == 'interval':
            raise ValueError("Table must be an interval table")
        super().set_table(table_name, table_data)

    def slice(self, slicer_table, target_table, target_se=None, columns=None, interval_include='both', sample_rate=None):
        """
        Slice the MSessionExplorer.Data object.
        """
        # Get the slicer table
        if isinstance(slicer_table, str):
            slicer_table = self.get_table(slicer_table)
        elif not isinstance(slicer_table, pd.DataFrame):
            raise ValueError("Slicer table must be a string or a pandas DataFrame")

        # Get the target table
        if isinstance(target_table, str):
            if target_se is None:
                target_se = self.target_se
            assert target_se is not None
            target_table = target_se.get_table(target_table)
        elif not isinstance(target_table, pd.DataFrame):
            raise ValueError("Target table must be a string or a pandas DataFrame")

        # Get the target table type
        target_table_type = self._get_table_type(target_table)

        # Select a subset of columns
        if columns is not None:
            # Make sure columns is a list
            if isinstance(columns, np.ndarray):
                columns = columns.tolist()
            elif not isinstance(columns, (tuple, list)):
                columns = [columns]
            
            # Convert column names to column indices
            columns = [x if isinstance(x, int) else target_table.columns.get_loc(x) for x in columns]

            # Must include time column for timeseries tables
            if target_table_type == 'timeseries':
                if 0 in columns:
                    columns.remove(0)
                columns = [0] + columns
            
            target_table = target_table.iloc[:, columns]

        # Slice the target table
        slices = []
        for t1, t2 in zip(slicer_table['start_time'], slicer_table['stop_time']):
            if target_table_type == 'interval':
                slices.append(self._slice_intervals(target_table, t1, t2, interval_include, sample_rate))

            elif target_table_type == 'timeseries':
                slices.append(self._slice_timeseries(target_table, t1, t2, sample_rate))

            elif target_table_type == 'timestamp':
                slices.append(self._slice_timestamps(target_table, t1, t2, sample_rate))

        # Concatenate the slices and reset indices
        slices = pd.concat(slices).reset_index(drop=True)

        return slices

    def _slice_timeseries(self, target_df, start_time, stop_time, sample_rate):
        """
        Slice a timeseries table.
        """
        # Masking
        mask = (target_df['time'] >= start_time) & (target_df['time'] <= stop_time)
        masked_df = target_df.loc[mask]

        # Make the masked table a single row
        row_df = {}
        for col in masked_df.columns:
            row_df[col] = [masked_df[col].values]
        row_df = pd.DataFrame(row_df)
        
        if sample_rate is None:
            return row_df
        
        # Make timestamps (avoid precision issues)
        t_query = np.arange(0, np.round(stop_time-start_time, 6), 1.0/sample_rate)
        t_query += start_time

        # Resample the sliced table to the sample rate
        resampled_df = {'time': [t_query]}
        for col in row_df.columns[1:]:
            resampled_df[col] = [np.interp(t_query, row_df['time'].iloc[0], row_df[col].iloc[0])]
        resampled_df = pd.DataFrame(resampled_df)

        return resampled_df

    def _slice_intervals(self, target_table, start_time, stop_time, interval_include, sample_rate):
        """
        Slice an interval table.
        """
        # Masking
        if interval_include == 'both':
            mask = (target_table['start_time'] >= start_time) & (target_table['stop_time'] <= stop_time)
        elif interval_include == 'start':
            mask = (target_table['start_time'] >= start_time) & (target_table['start_time'] <= stop_time)
        elif interval_include == 'stop':
            mask = (target_table['stop_time'] >= start_time) & (target_table['stop_time'] <= stop_time)
        elif interval_include == 'either':
            mask = (target_table['stop_time'] >= start_time) & (target_table['start_time'] <= stop_time)
        else:
            raise ValueError(f"Invalid interval_include value: {interval_include}")
        masked_table = target_table.loc[mask]
        
        if sample_rate is None:
            # Make the masked table a single row
            row_table = {}
            for col in masked_table.columns:
                row_table[col] = [masked_table[col].values]
            row_table = pd.DataFrame(row_table)
            return row_table
        else:
            # Construct boxcar timeseries table
            timestamps = np.arange(start_time, stop_time, 1.0/sample_rate)
            boxcar = np.zeros(len(timestamps))
            for _, row in masked_table.iterrows():
                mask = (timestamps >= row['start_time']) & (timestamps < row['stop_time'])
                boxcar[mask] += 1
            boxcar_table = pd.DataFrame({'time': [timestamps], 'boxcar': [boxcar]})
            return boxcar_table
    
    def _slice_timestamps(self, target_table, start_time, stop_time, sample_rate):
        """
        Slice a timestamp column.
        """
        if sample_rate is None:
            output_table = {}
        else:
            # Create output table with timestamps at desired sample rate
            timestamps = np.arange(start_time, stop_time, 1.0/sample_rate)
            output_table = {'time': [timestamps]}

        # Process each column of timestamps
        for col in target_table.columns:
            # Get timestamps for this column
            col_timestamps = target_table[col].iloc[0]
            
            # Mask timestamps to time window
            mask = (col_timestamps >= start_time) & (col_timestamps < stop_time)
            masked_timestamps = col_timestamps[mask]

            if sample_rate is None:
                # Just return the masked timestamps
                output_table[col] = [masked_timestamps]
            else:
                # Create impulse timeseries using histogram
                impulse, _ = np.histogram(masked_timestamps, bins=timestamps)
                impulse[-2] += impulse[-1]  # add the last bin's count to the second-to-last bin since histogram excludes the right edge
                impulse = impulse[:-1]
                output_table[col] = [impulse]

        return pd.DataFrame(output_table)

    def align_time(self, table, reference_times):
        """
        Align the time of a table to new reference times.
        """
        if isinstance(reference_times, pd.Series):
            reference_times = reference_times.values
        if len(reference_times) == 1:
            reference_times = [reference_times] * len(table)
        elif len(reference_times) != len(table):
            raise ValueError("Reference times must be a scalar or a list/array of the same length as the table")

        table_type = self._get_table_type(table)
        if table_type == 'timeseries':
            table['time'] = [t - t0 for t, t0 in zip(table['time'].values, reference_times)]
        elif table_type == 'interval':
            table['start_time'] = [t - t0 for t, t0 in zip(table['start_time'].values, reference_times)]
            table['stop_time'] = [t - t0 for t, t0 in zip(table['stop_time'].values, reference_times)]
        elif table_type == 'timestamp':
            for col in table.columns:
                table[col] = [t - t0 for t, t0 in zip(table[col].values, reference_times)]
        
        return table
