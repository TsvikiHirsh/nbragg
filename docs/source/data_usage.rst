.. data_usage:

Data Handling in nres
=====================

The nres package utilizes pandas DataFrame as its core data structure for handling neutron transmission data. This approach leverages the powerful data manipulation and analysis capabilities of pandas, providing flexibility and efficiency in working with experimental data.

Data Class
----------

The primary interface for data handling in nres is the `Data` class. This class encapsulates methods for reading, processing, and visualizing neutron transmission data.

Key Attributes:
^^^^^^^^^^^^^^^

- ``table``: A pandas DataFrame containing energy, transmission, and error values.
- ``tgrid``: A pandas Series representing the time-of-flight grid corresponding to the data.

Data Import Methods
-------------------

The `Data` class provides two main methods for importing data:

1. ``from_counts``:
   This method creates a `Data` object from raw counts data, typically from signal and open beam measurements.

   .. code-block:: python

      data = Data.from_counts(signal="signal.csv", openbeam="openbeam.csv")

   This method:
   - Reads counts data from CSV files
   - Calculates transmission and associated errors
   - Converts time-of-flight to energy
   - Applies background correction if empty signal/openbeam data is provided

2. ``from_transmission``:
   This method creates a `Data` object directly from a transmission data file.

   .. code-block:: python

      data = Data.from_transmission("transmission_data.txt")

   This method reads a file containing energy, transmission, and error values, typically space-separated.

Pandas Integration
------------------

By using pandas DataFrame as the underlying data structure, nres allows users to leverage the full power of pandas for data manipulation and analysis. Some key advantages include:

1. **Data Filtering**: 
   You can easily filter your data based on energy ranges or other criteria:

   .. code-block:: python

      # Filter data for energies between 1 eV and 1000 eV
      filtered_data = data.table.query("1 <= energy <= 1000")

2. **Uncertainty Handling**:
   Manage uncertainty values with pandas operations:

   .. code-block:: python

      # Remove data points with relative uncertainty > 10%
      clean_data = data.table[data.table['err'] / data.table['trans'] <= 0.1]

3. **Data Transformation**:
   Apply complex transformations to your data:

   .. code-block:: python

      # Add a column for relative error
      data.table['rel_err'] = data.table['err'] / data.table['trans']

4. **Statistical Analysis**:
   Utilize pandas' statistical functions:

   .. code-block:: python

      # Calculate mean transmission in specific energy range
      mean_trans = data.table.query("1e6 <= energy <= 2e6")['trans'].mean()

5. **Data Export**:
   Easily export your data to various formats:

   .. code-block:: python

      # Export to CSV
      data.table.to_csv("processed_data.csv", index=False)

Visualization
-------------

The `Data` class includes a `plot` method for quick visualization of the transmission data:

.. code-block:: python

   data.plot(xlim=(0.5e6, 1e7), ylim=(0, 1), logx=True)

This method utilizes pandas' plotting capabilities, which are built on matplotlib, allowing for easy customization of plots.

Advanced Usage
--------------

The integration with pandas allows for advanced data handling techniques:

1. **Merging Datasets**:
   Combine data from multiple experiments:

   .. code-block:: python

      merged_data = pd.concat([data1.table, data2.table]).sort_values('energy')

2. **Rolling Statistics**:
   Apply rolling window calculations:

   .. code-block:: python

      data.table['rolling_mean'] = data.table['trans'].rolling(window=10).mean()

3. **Grouping and Aggregation**:
   Group data by energy bins and perform aggregations:

   .. code-block:: python

      binned_data = data.table.groupby(pd.cut(data.table['energy'], bins=100)).agg({
          'trans': 'mean',
          'err': lambda x: np.sqrt(np.sum(x**2)) / len(x)
      })

Conclusion
----------

The use of pandas DataFrame in nres provides a powerful and flexible foundation for handling neutron transmission data. It allows users to perform complex data operations, statistical analysis, and visualizations with ease, while maintaining the specific requirements of neutron resonance spectroscopy data processing.

By leveraging pandas, nres combines the specificity needed for neutron data analysis with the broad capabilities of a leading data manipulation library, offering a robust toolset for researchers in the field.