import pandas as pd
from fuzzywuzzy import process
import warnings
import numpy as np
import pandas as pd
from typing import List, Dict, Tuple, Optional, Union, Any, Sequence
from joblib import Parallel, delayed
import statsmodels.api as sm
warnings.filterwarnings("ignore", category=FutureWarning)


################ functions for search the valid medications ################
def fuzzy_search(cell, word_list, threshold=90):
    if pd.isna(cell):
        return False
    cell_str = str(cell).strip()
    if not cell_str or cell_str in {'?', '-', '--', 'N/A', 'na'}:
        return False
    for word in word_list:
        match = process.extractOne(cell_str, [word])
        if match and match[1] >= threshold:
            return True
    return False

def contains_any_drug(text, drug_list):
    return any(drug.lower() in str(text).lower() for drug in drug_list)


def get_matched_values(df, drug_list, threshold=90, n_jobs=-1):
    # flatten df into list of (row, col, value)
    flat = [(i, j, df.iat[i, j]) for i in range(df.shape[0]) for j in range(df.shape[1])]

    # run fuzzy_search in parallel
    results = Parallel(n_jobs=n_jobs)(
        delayed(fuzzy_search)(val, drug_list, threshold) for _, _, val in flat
    )

    # collect matched values
    matched_values = [val for (match, (_, _, val)) in zip(results, flat) if match]
    return list(set(matched_values))



def extract_medication_info(
    df: pd.DataFrame,
    drug_values: Sequence[str],
    medication_name_columns: Sequence[str],
    selected_columns: Sequence[str],
    *,
    id_col: str = "id",
    visit_col: str = "visit_no",
) -> pd.DataFrame:
    """
    Keep visits where ANY medication name matches one of `drug_values`.
    Returns first matching Medication_Name per visit + id/visit + selected columns.
    """
    # Normalize lookup set (trim + lower)
    drug_set = {str(x).strip().lower() for x in drug_values}

    # Use only med columns that exist
    med_cols = [c for c in medication_name_columns if c in df.columns]
    if not med_cols:
        return pd.DataFrame(columns=["Medication_Name", id_col, visit_col, *selected_columns])

    # Row-wise: does any med column equal a value in `drug_set`?
    def has_match_row(row) -> bool:
        for c in med_cols:
            v = row.get(c)
            if pd.isna(v):
                continue
            if str(v).strip().lower() in drug_set:
                return True
        return False

    filtered_df = df[df.apply(has_match_row, axis=1)].copy()
    if filtered_df.empty:
        return pd.DataFrame(columns=["Medication_Name", id_col, visit_col, *selected_columns])

    # First matching medication name per row
    def first_match_name(row):
        for c in med_cols:
            v = row.get(c)
            if pd.isna(v):
                continue
            if str(v).strip().lower() in drug_set:
                return v
        return None

    med_name = filtered_df.apply(first_match_name, axis=1).rename("Medication_Name").to_frame()

    # Carry id/visit + selected columns (dedup + existing only)
    carry_cols = [id_col, visit_col, *selected_columns]
    carry_cols = list(dict.fromkeys(carry_cols))  # de-dup
    carry_cols = [c for c in carry_cols if c in filtered_df.columns]

    out = pd.concat([med_name.reset_index(drop=True),
                     filtered_df[carry_cols].reset_index(drop=True)], axis=1)
    return out.reset_index(drop=True)


def medication_all_timepoints(
    df: pd.DataFrame,
    drug_values: Sequence[str],
    medication_name_columns: Sequence[str],
    selected_columns: Sequence[str],
    *,
    id_col: str = "id",
    visit_col: str = "visit_no",
) -> pd.DataFrame:
    """
    Return all visits for participants who EVER have a target medication, with:
      - Has_Medication_This_Visit (bool)
      - Has_Medication_Ever (cumulative within id)
    Also adds first matching Medication_Name for each visit.
    """
    # Normalize lookup set (trim + lower)
    drug_set = {str(x).strip().lower() for x in drug_values}

    # Identify visits with any match
    med_cols = [c for c in medication_name_columns if c in df.columns]
    if not med_cols:
        base_cols = [id_col, visit_col, *selected_columns]
        out = df.head(0)[[c for c in base_cols if c in df.columns]].copy()
        out["Has_Medication_This_Visit"] = pd.Series(dtype=bool)
        out["Has_Medication_Ever"] = pd.Series(dtype=bool)
        out["Medication_Name"] = pd.Series(dtype=object)
        return out

    def row_has_match(row) -> bool:
        for c in med_cols:
            v = row.get(c)
            if pd.isna(v):
                continue
            if str(v).strip().lower() in drug_set:
                return True
        return False

    def first_match_name(row):
        for c in med_cols:
            v = row.get(c)
            if pd.isna(v):
                continue
            if str(v).strip().lower() in drug_set:
                return v
        return None

    # People who ever match
    ever_mask = df.apply(row_has_match, axis=1)
    ids_with_med = df.loc[ever_mask, id_col].unique()
    df_filtered = df[df[id_col].isin(ids_with_med)].copy()

    # Flags per visit
    df_filtered["Has_Medication_This_Visit"] = df_filtered.apply(row_has_match, axis=1)

    # Sort by visit for cumulative ever flag
    # visit_sort = pd.to_numeric(df_filtered[visit_col], errors="coerce")
    # df_filtered = df_filtered.assign(_vsort=visit_sort)
    # df_filtered = df_filtered.sort_values([id_col, "_vsort", visit_col]).drop(columns="_vsort")
    df_filtered = df_filtered.sort_values([id_col, visit_col])

    df_filtered["Has_Medication_Ever"] = df_filtered.groupby(id_col)["Has_Medication_This_Visit"].cummax()

    # Per-visit first matching med name
    df_filtered["Medication_Name"] = df_filtered.apply(first_match_name, axis=1)
    df_filtered = df_filtered.drop(columns=[c for c in medication_name_columns if c in df_filtered.columns])
    return df_filtered.reset_index(drop=True)



########### for control group
def contains_any_medication_type(cell, drug_set):
    if pd.isna(cell):
        return False
    text = str(cell).lower()
    return any(drug.lower() in text for drug in drug_set)


def rename_and_select_columns(df, excel_path, old_col_name, new_col_name, extra_cols=None):
    mapping_df = pd.read_excel(excel_path)
    mapping_dict = dict(zip(mapping_df[old_col_name], mapping_df[new_col_name]))
    df = df.rename(columns=mapping_dict)
    keep_cols = list(mapping_dict.values()) + (extra_cols or [])
    return df[keep_cols]


def merge_medication_longitudinal(datasets: dict, control_key: str | None = None, extra_cols: list[str] | None = None) -> pd.DataFrame:

    base_cols = ['id','visit_date','status','age','sex','edu','APOE4']
    med_keys = [k for k in ['ACEi','ARB','BetaBlk','CCB','Diuretic','Statin','Metformin'] if k in datasets]

    # 1) Master visit grid
    master = (pd.concat([d[base_cols + (extra_cols or [])] for d in datasets.values()], ignore_index=True)
                .assign(visit_date=lambda x: pd.to_datetime(x['visit_date'], errors='coerce')))
    master = (master.groupby(['id','visit_date'], as_index=False)
                    .agg({c:(lambda s: s.dropna().iloc[0] if s.notna().any() else np.nan) 
                          for c in base_cols + (extra_cols or []) if c not in ['id','visit_date']}))

    # 2) Add medication flags
    for mk in med_keys:
        df = datasets[mk].copy()
        df['visit_date'] = pd.to_datetime(df['visit_date'], errors='coerce')
        # pos = df[df['Medication'].fillna(False).astype(bool)][['id','visit_date']].drop_duplicates()
        pos = df[df['Has_Medication_This_Visit'].fillna(False).astype(bool)][['id','visit_date']].drop_duplicates()
        master = master.merge(pos.assign(**{mk: True}), on=['id','visit_date'], how='left')
        master[mk] = master[mk].fillna(False)

    # 3) Control rows
    if control_key and control_key in datasets:
        ctrl = datasets[control_key][['id','visit_date']].copy()
        ctrl['visit_date'] = pd.to_datetime(ctrl['visit_date'], errors='coerce')
        master = master.merge(ctrl.assign(__ctrl=1), on=['id','visit_date'], how='left')
        for mk in med_keys:
            master.loc[master['__ctrl'].notna(), mk] = False
        master = master.drop(columns='__ctrl')

    # 4) Final tidy columns
    # master = master.rename(columns={'edu':'education','APOE':'APOE4'})
    master['Total_Meds'] = master[med_keys].astype(int).sum(axis=1) if med_keys else 0

    core_order = ['id','visit_date','status','age','sex','edu','APOE4'] + med_keys + ['Total_Meds']
    if extra_cols:
        core_order += [c for c in extra_cols if c in master.columns]

    return master[core_order].sort_values(['id','visit_date']).reset_index(drop=True)



################################ HABSHD ################################
def clean_dataframe_advanced(df, missing_values=None, threshold=0.7, 
                            drop_rows_threshold=None, verbose=True):
    
    if missing_values is None:
        missing_values = [-777777, -777777.0, -9999.000000, -9999, -999, -999.0, -8888, "-9999"]
    
    df_cleaned = df.copy()
    
    # Replace missing values
    df_cleaned = df_cleaned.replace(missing_values, np.nan)
    
    # Drop columns with high missing rates
    col_missing_rates = df_cleaned.isnull().sum() / len(df_cleaned)
    cols_to_drop = col_missing_rates[col_missing_rates > threshold].index.tolist()
    df_cleaned = df_cleaned.drop(columns=cols_to_drop)
    
    # Optionally drop rows with high missing rates
    rows_dropped = 0
    if drop_rows_threshold is not None:
        before_rows = len(df_cleaned)
        row_missing_rates = df_cleaned.isnull().sum(axis=1) / len(df_cleaned.columns)
        df_cleaned = df_cleaned[row_missing_rates <= drop_rows_threshold]
        rows_dropped = before_rows - len(df_cleaned)
    
    if verbose:
        print(f"Original shape: {df.shape}")
        print(f"Missing values replaced: {missing_values}")
        print(f"Columns dropped: {len(cols_to_drop)} ({cols_to_drop})")
        if drop_rows_threshold is not None:
            print(f"Rows dropped: {rows_dropped}")
        print(f"Final shape: {df_cleaned.shape}")
    return df_cleaned

def merge_dataframes(dataframes: List[pd.DataFrame], 
                    keys: List[str] = ['Visit_ID', 'Med_ID'],
                    how: str = 'outer',
                    handle_duplicates: str = 'auto',
                    df_names: Optional[List[str]] = None,
                    verbose: bool = True) -> pd.DataFrame:
    
    if len(dataframes) < 2:
        raise ValueError("At least 2 DataFrames are required for merging")
    
    if df_names is None:
        df_names = [f"df_{i}" for i in range(len(dataframes))]
    elif len(df_names) != len(dataframes):
        raise ValueError("Number of df_names must match number of dataframes")
    
    # Validate merge keys exist in all DataFrames
    for i, df in enumerate(dataframes):
        missing_keys = [key for key in keys if key not in df.columns]
        if missing_keys:
            raise ValueError(f"DataFrame {i} ({df_names[i]}) missing keys: {missing_keys}")
    
    # Detect column conflicts
    conflicts = detect_column_conflicts(dataframes, keys, df_names)
    
    if verbose and conflicts:
        print("Column conflicts detected:")
        for col, dfs in conflicts.items():
            print(f"  '{col}' appears in: {dfs}")
    
    processed_dfs = handle_column_conflicts(dataframes, conflicts, handle_duplicates, df_names, keys)
    result = processed_dfs[0]
    merge_info = []
    
    for i, df in enumerate(processed_dfs[1:], 1):
        before_rows = len(result)
        result = result.merge(df, on=keys, how=how, suffixes=('', f'_from_{df_names[i]}'))
        after_rows = len(result)
        merge_info.append({
            'step': i,
            'df_name': df_names[i],
            'rows_before': before_rows,
            'rows_after': after_rows,
            'columns_added': len(df.columns) - len(keys)
        })
    
    if verbose:
        print(f"\nMerge completed:")
        print(f"Final DataFrame shape: {result.shape}")
        print(f"Merge key(s): {keys}")
        print(f"Merge type: {how}")
        
        print("\nMerge steps:")
        for info in merge_info:
            print(f"  Step {info['step']}: Added {df_names[info['step']]} "
                  f"({info['rows_before']} -> {info['rows_after']} rows, "
                  f"+{info['columns_added']} columns)")
    return result

def detect_column_conflicts(dataframes: List[pd.DataFrame], 
                          keys: List[str], 
                          df_names: List[str]) -> Dict[str, List[str]]:
    conflicts = {}
    for i, df in enumerate(dataframes):
        non_key_cols = [col for col in df.columns if col not in keys]
        for col in non_key_cols:
            if col not in conflicts:
                conflicts[col] = []
            conflicts[col].append(df_names[i])
    conflicts = {col: dfs for col, dfs in conflicts.items() if len(dfs) > 1}
    return conflicts

def handle_column_conflicts(dataframes: List[pd.DataFrame], 
                          conflicts: Dict[str, List[str]], 
                          strategy: str,
                          df_names: List[str],
                          keys: List[str]) -> List[pd.DataFrame]:
    processed_dfs = [df.copy() for df in dataframes]
    if strategy == 'auto' or strategy == 'compare':
        if conflicts:
            print(f"\nAnalyzing {len(conflicts)} column conflicts...")
            conflict_analysis = analyze_conflicts(dataframes, conflicts, df_names, keys)
            
            if strategy == 'compare':
                print_conflict_comparison(conflict_analysis)
                return processed_dfs
            processed_dfs = auto_resolve_conflicts(processed_dfs, conflict_analysis, df_names)
    
    elif strategy == 'rename':
        for i, df in enumerate(processed_dfs):
            for col in conflicts:
                if col in df.columns:
                    new_name = f"{col}_{df_names[i]}"
                    df.rename(columns={col: new_name}, inplace=True)
    
    elif strategy == 'keep_first':
        for i, df in enumerate(processed_dfs[1:], 1):
            cols_to_drop = [col for col in conflicts if col in df.columns]
            df.drop(columns=cols_to_drop, inplace=True)
    
    elif strategy == 'keep_last':
        for i, df in enumerate(processed_dfs[:-1]):
            cols_to_drop = [col for col in conflicts if col in df.columns]
            df.drop(columns=cols_to_drop, inplace=True)
    return processed_dfs

def analyze_conflicts(dataframes: List[pd.DataFrame], 
                     conflicts: Dict[str, List[str]], 
                     df_names: List[str],
                     keys: List[str]) -> Dict:
    analysis = {}
    
    for col in conflicts:
        col_analysis = {
            'column': col,
            'dataframes': conflicts[col],
            'identical_values': True,
            'data_types': [],
            'null_counts': [],
            'unique_values': [],
            'recommendation': 'keep_any'
        }
        
        # Get data from each DataFrame for this column
        col_data = []
        for i, df in enumerate(dataframes):
            if col in df.columns:
                # Merge with keys to align data
                temp_df = df[keys + [col]].copy()
                col_data.append((df_names[i], temp_df))
                col_analysis['data_types'].append(str(df[col].dtype))
                col_analysis['null_counts'].append(df[col].isnull().sum())
                col_analysis['unique_values'].append(df[col].nunique())
        
        # Check if values are identical across DataFrames
        if len(col_data) > 1:
            # Merge all DataFrames on keys to compare values
            merged_for_comparison = col_data[0][1]
            for name, data in col_data[1:]:
                merged_for_comparison = merged_for_comparison.merge(
                    data, on=keys, how='outer', suffixes=('', f'_{name}')
                )
            
            # Check for differences
            comparison_cols = [c for c in merged_for_comparison.columns if c not in keys]
            if len(comparison_cols) > 1:
                # Compare values
                base_col = comparison_cols[0]
                for other_col in comparison_cols[1:]:
                    mask = merged_for_comparison[base_col].notna() & merged_for_comparison[other_col].notna()
                    if mask.any():
                        differences = merged_for_comparison[mask][base_col] != merged_for_comparison[mask][other_col]
                        if differences.any():
                            col_analysis['identical_values'] = False
                            col_analysis['recommendation'] = 'rename_all'
                            break
        
        analysis[col] = col_analysis
    
    return analysis

def auto_resolve_conflicts(dataframes: List[pd.DataFrame], 
                          analysis: Dict, 
                          df_names: List[str]) -> List[pd.DataFrame]:
    processed_dfs = [df.copy() for df in dataframes]
    
    for col, info in analysis.items():
        if info['recommendation'] == 'rename_all':
            # Rename in all DataFrames
            for i, df in enumerate(processed_dfs):
                if col in df.columns:
                    new_name = f"{col}_{df_names[i]}"
                    df.rename(columns={col: new_name}, inplace=True)
            print(f"  Renamed '{col}' in all DataFrames (values differ)")
        
        elif info['recommendation'] == 'keep_any':
            # Keep from first DataFrame, drop from others
            for i, df in enumerate(processed_dfs[1:], 1):
                if col in df.columns:
                    df.drop(columns=[col], inplace=True)
            print(f"  Kept '{col}' from first DataFrame (values identical)")
    
    return processed_dfs

def print_conflict_comparison(analysis: Dict):
    """
    Print detailed comparison of conflicting columns.
    """
    print("\nDetailed Conflict Analysis:")
    print("=" * 60)
    
    for col, info in analysis.items():
        print(f"\nColumn: '{col}'")
        print(f"  Appears in: {', '.join(info['dataframes'])}")
        print(f"  Data types: {info['data_types']}")
        print(f"  Null counts: {info['null_counts']}")
        print(f"  Unique values: {info['unique_values']}")
        print(f"  Identical values: {info['identical_values']}")
        print(f"  Recommendation: {info['recommendation']}")

################# NACC #################
def add_age_column_NACC(df):
    # Use existing NACC_VISDATE column (ensure it's datetime)
    df['VISITDATE'] = pd.to_datetime(df['NACC_VISDATE'], errors='coerce')
    
    # Construct birth date from birth year and month (only for valid numeric values)
    valid_birth_mask = (
        pd.notna(df['BIRTHYR']) & 
        pd.notna(df['BIRTHMO']) & 
        (df['BIRTHYR'] > 0) & 
        (df['BIRTHMO'] > 0)
    )
    
    df['BIRTHDATE'] = pd.NaT
    df.loc[valid_birth_mask, 'BIRTHDATE'] = pd.to_datetime(
        dict(
            year=df.loc[valid_birth_mask, 'BIRTHYR'], 
            month=df.loc[valid_birth_mask, 'BIRTHMO'], 
            day=1
        ),
        errors='coerce'
    )
    
    # Calculate age in years (only where both dates are valid)
    valid_date_mask = pd.notna(df['VISITDATE']) & pd.notna(df['BIRTHDATE'])
    df['AGE'] = np.nan  # Use np.nan instead of pd.NA
    df.loc[valid_date_mask, 'AGE'] = (
        (df.loc[valid_date_mask, 'VISITDATE'] - df.loc[valid_date_mask, 'BIRTHDATE']).dt.total_seconds() 
        / (365.25 * 24 * 60 * 60)
    ).round(2)
    
    # Debug: Print data validity info
    print(f"Valid visit dates: {pd.notna(df['VISITDATE']).sum()}")
    print(f"Valid birth dates: {pd.notna(df['BIRTHDATE']).sum()}")
    print(f"Valid calculated ages: {pd.notna(df['AGE']).sum()}")
    print(f"Sample AGE values: {df['AGE'].head(10).tolist()}")
    
    # Process NACCAGE if available
    if 'NACCAGE' in df.columns:
        print("NACCAGE column found - processing age comparisons and imputations")
        
        # Count missing NACCAGE values before imputation
        missing_naccage_before = df['NACCAGE'].isna().sum()
        
        # Fill missing NACCAGE with calculated AGE
        df['NACCAGE'].fillna(df['AGE'], inplace=True)
        
        # Count how many NACCAGE were imputed with AGE
        missing_naccage_after = df['NACCAGE'].isna().sum()
        imputed_count = missing_naccage_before - missing_naccage_after
        print(f"Missing NACCAGE values imputed with calculated AGE: {imputed_count}")
        
        # Replace AGE with NACCAGE where they differ by less than 1 year and both are valid
        valid_both_mask = pd.notna(df['AGE']) & pd.notna(df['NACCAGE'])
        close_match_mask = valid_both_mask & (abs(df['AGE'] - df['NACCAGE']) < 1)
        
        replaced_count = close_match_mask.sum()
        df.loc[close_match_mask, 'NACCAGE'] = df.loc[close_match_mask, 'AGE']
        print(f"AGE values replaced with NACCAGE (difference < 1 year): {replaced_count}")
    else:
        print("NACCAGE column not found - using only calculated AGE")
    
    return df



def clean_status_and_demos(df, id_col="id", visit_col="visit_no", status_col="status", fill_cols=None):
    if fill_cols is None: 
        fill_cols = ["sex","edu","APOE4"]

    df = df.sort_values([id_col, visit_col]).copy()
    s = df[status_col].copy(); missing_before = s.isna().sum()
    hc_fill = s.where(s.eq("HC")).groupby(df[id_col]).bfill()
    s = s.where(~s.isna(), hc_fill)
    ad_seen = s.eq("AD").groupby(df[id_col]).cumsum().gt(0)
    df[status_col] = np.where(ad_seen, "AD", s)

    missing_after = df[status_col].isna().sum()
    cols_to_fill = [c for c in (fill_cols or []) if c in df.columns]
    if cols_to_fill:
        df[cols_to_fill] = df.groupby(id_col)[cols_to_fill].transform(lambda x: x.ffill().bfill())

    print(f"Filled {missing_before - missing_after} missing {status_col} records; columns forward/back filled: {cols_to_fill}")
    return df

#################### AIBL ####################
def fill_missing_ages_AIBL(df, interval_years=1.5):
    """
    Fill missing age values based on 18-month intervals between collections.
    
    Parameters:
    df: DataFrame with columns 'AIBL ID', 'Collection', 'Age'
    interval_years: Time interval between collections in years (default 1.5)
    
    Returns:
    DataFrame with missing ages filled
    """
    
    # Create a copy to avoid modifying the original dataframe
    df_filled = df.copy()
    
    # Sort by AIBL ID and Collection to ensure proper order
    df_filled = df_filled.sort_values(['AIBL ID', 'Collection']).reset_index(drop=True)
    
    # Group by AIBL ID to process each individual separately
    def fill_individual_ages(group):
        group = group.copy()
        group = group.sort_values('Collection').reset_index(drop=True)
        
        # Find non-missing ages to use as reference points
        valid_ages = group.dropna(subset=['Age'])
        
        if len(valid_ages) == 0:
            # No valid ages for this individual
            return group
        
        # For each missing age, find the nearest known age and calculate
        for idx in group.index:
            if pd.isna(group.loc[idx, 'Age']):
                current_collection = group.loc[idx, 'Collection']
                
                # Try forward fill first (from earlier collections)
                earlier_ages = valid_ages[valid_ages['Collection'] < current_collection]
                if len(earlier_ages) > 0:
                    # Use the most recent earlier age
                    ref_row = earlier_ages.iloc[-1]
                    collection_diff = current_collection - ref_row['Collection']
                    estimated_age = ref_row['Age'] + collection_diff * interval_years
                    group.loc[idx, 'Age'] = estimated_age
                    continue
                
                # Try backward fill (from later collections)
                later_ages = valid_ages[valid_ages['Collection'] > current_collection]
                if len(later_ages) > 0:
                    # Use the earliest later age
                    ref_row = later_ages.iloc[0]
                    collection_diff = ref_row['Collection'] - current_collection
                    estimated_age = ref_row['Age'] - collection_diff * interval_years
                    group.loc[idx, 'Age'] = estimated_age
        
        return group
    
    # Apply the filling function to each individual
    df_filled = df_filled.groupby('AIBL ID').apply(fill_individual_ages).reset_index(drop=True)
    
    return df_filled


############# for preparing the dataset in R 
def clean_and_filter_participants(
    df: pd.DataFrame,
    *,
    id_col: str = 'AIBL ID',
    class_col: str = 'Neuropsych.Simple Classification',
    medication_ever_col: Optional[str] = "Has_Medication_Ever",
    time_col: str = 'Collection',
    allowed_classes = ('HC', 'MCI', 'AD', 'MCIX'),
    min_visits: int = 2,
    if_print: bool = True
) -> pd.DataFrame:
    """
    1.  Keep only rows whose `class_col` is in allowed_classes.
        • Convert 'MCIX' → 'MCI'.
    2.  Keep participants who have **≥ min_visits** rows remaining.
    3.  Keep participants who have **at least one** row with
        `medication_ever_col == 1`.
    4.  Optionally print how many participants remain & were dropped.
    
    Returns
    -------
    cleaned_df : pd.DataFrame
        The filtered dataframe.
    """
    start_ids = set(df[id_col].unique())
    df1 = df[df[class_col].isin(allowed_classes)].copy()
    df1.loc[df1[class_col] == 'MCIX', class_col] = 'MCI'

    # Drop participants who have HC or MCI after being diagnosed with AD
    ids_to_drop = []
    for pid, grp in df1.groupby(id_col):
        ad_times = grp.loc[grp[class_col] == 'AD', time_col]
        if not ad_times.empty:
            first_ad_time = ad_times.min()
            if any((grp[time_col] > first_ad_time) & (grp[class_col].isin(['HC', 'MCI']))):
                ids_to_drop.append(pid)
    df1 = df1[~df1[id_col].isin(ids_to_drop)]

    visit_counts = df1.groupby(id_col)[class_col].size()
    ids_enough_visits = visit_counts[visit_counts >= min_visits].index
    df2 = df1[df1[id_col].isin(ids_enough_visits)].copy()


    # ids_with_med = (
    #     df2.loc[df2[medication_ever_col] == 1, id_col]
    #         .unique()
    # )
    # cleaned_df = df2[df2[id_col].isin(ids_with_med)].copy()
    if medication_ever_col and medication_ever_col in df2.columns:
        # coerce to numeric if needed
        if not pd.api.types.is_numeric_dtype(df2[medication_ever_col]):
            df2[medication_ever_col] = pd.to_numeric(
                df2[medication_ever_col], errors="coerce"
            )

        ids_with_med = df2.loc[
            df2[medication_ever_col] == 1, id_col
        ].unique()
        cleaned_df = df2[df2[id_col].isin(ids_with_med)].copy()
    else:
        # no medication filter applied
        cleaned_df = df2.copy()

    if if_print:
        end_ids   = set(cleaned_df[id_col].unique())
        dropped   = start_ids - end_ids
        print(f'Remaining participants: {len(end_ids):,}')
        print(f'Dropped participants:   {len(dropped):,}')
    return cleaned_df


def clean_longitudinal(df, id_col='id', time_col='months_since_baseline', visit_col='visit_no'):
    df = df.copy()
    # Sort by id and time
    df.sort_values([id_col, time_col], kind='mergesort', inplace=True)
    # Re-base months_since_baseline so earliest time is 0
    t0 = df.groupby(id_col)[time_col].transform('min')
    df[time_col] = df[time_col] - t0
    # Re-number visit_no starting from 1
    df[visit_col] = df.groupby(id_col)[time_col].rank(method='first').astype(int)
    return df.reset_index(drop=True)

# def baseline_summary(df):
#     df = df.copy()
#     # pick baseline rows (visit_no==1 preferred, else months_since_baseline==0)
#     baseline = df[df['visit_no'] == 1] if 'visit_no' in df else df[df['months_since_baseline'] == 0]
#     baseline = baseline.sort_values(['id', 'months_since_baseline']).drop_duplicates('id', keep='first')

#     n = baseline['id'].nunique()  # number of participants

#     # helper formatters
#     fmt_mean_sd = lambda s: f"{s.mean():.1f} ± {s.std(ddof=1):.1f}" if len(s)>0 else np.nan
#     fmt_count_pct = lambda s: f"{s} ({s/n*100:.1f}%)" if n>0 else "0 (0.0%)"
#     fmt_iqr = lambda s: f"{s.median():.1f} [{s.quantile(0.25):.1f}–{s.quantile(0.75):.1f}]" if len(s)>0 else np.nan
#     fmt_minmax = lambda s: f"{s.min():.1f} – {s.max():.1f}" if len(s)>0 else np.nan

#     # metrics
#     age = fmt_mean_sd(baseline['age'])
#     female = fmt_count_pct(((baseline['sex'].astype(str).str.upper().isin(['F','FEMALE'])) | 
#                         (baseline['sex'] == 1)).sum())
#     edu = fmt_mean_sd(baseline['edu'])
#     # apoe = fmt_count_pct((baseline['APOE4'].astype(str).str.upper().isin(['1','YES','TRUE',])).sum())
#     apoe = fmt_count_pct(((baseline['APOE4'].astype(str).str.upper().isin(['1','YES','TRUE'])) | 
#                         (baseline['APOE4'] == 1)).sum())
#     # visits = fmt_iqr(df.groupby('id').size())
#     visits = fmt_minmax(df.groupby('id').size())
#     # followup = fmt_minmax(df.groupby('id')['months_since_baseline'].max() - df.groupby('id')['months_since_baseline'].min())
#     followup = fmt_minmax(df.groupby('id')['months_since_baseline'].agg(lambda s: s.max(skipna=True) - s.min(skipna=True)).dropna())

#     cu = fmt_count_pct((baseline['status'].astype(str).str.upper() == 'HC').sum())
#     mci = fmt_count_pct((baseline['status'].astype(str).str.upper() == 'MCI').sum())
#     ad = fmt_count_pct((baseline['status'].astype(str).str.upper() == 'AD').sum())
#     meds = fmt_mean_sd(baseline['Total_Meds'])

#     cdr_change  = fmt_mean_sd(df.sort_values(['id','months_since_baseline']).groupby('id')['CDR'] .apply(lambda s: s.dropna().iloc[-1]-s.dropna().iloc[0] if s.dropna().size>1 else np.nan).dropna())
#     mmse_change = fmt_mean_sd(df.sort_values(['id','months_since_baseline']).groupby('id')['MMSE'].apply(lambda s: s.dropna().iloc[-1]-s.dropna().iloc[0] if s.dropna().size>1 else np.nan).dropna())

#     return pd.DataFrame({
#         "Measure": [
#             "Age at baseline (year)",
#             "Gender (Female)",
#             "Education (year)",
#             "APOE4 (YES)**",
#             "Visits (Record)",
#             "Follow up intervals (Month)",
#             "CU at baseline",
#             "MCI at baseline",
#             "AD at baseline",
#             "Average medication taken at baseline",
#             "Overall CDR change",
#             "Overall MMSE change"
#         ],
#         "Value": [
#             age, female, edu, apoe, visits, followup, cu, mci, ad, meds, cdr_change,mmse_change
#         ]
#     })
def baseline_summary(df):
    df = df.copy()
    # pick baseline rows (visit_no==1 preferred, else months_since_baseline==0)
    baseline = df[df['visit_no'] == 1] if 'visit_no' in df else df[df['months_since_baseline'] == 0]
    baseline = baseline.sort_values(['id', 'months_since_baseline']).drop_duplicates('id', keep='first')

    n = baseline['id'].nunique()  # number of participants

    # helper formatters
    fmt_mean_sd = lambda s: f"{s.mean():.1f} ± {s.std(ddof=1):.1f}" if len(s)>0 else np.nan
    fmt_count_pct = lambda s: f"{s} ({s/n*100:.1f}%)" if n>0 else "0 (0.0%)"
    fmt_iqr = lambda s: f"{s.median():.1f} [{s.quantile(0.25):.1f}–{s.quantile(0.75):.1f}]" if len(s)>0 else np.nan
    fmt_minmax = lambda s: f"{s.min():.1f} – {s.max():.1f}" if len(s)>0 else np.nan

    # metrics
    age = fmt_mean_sd(baseline['age'])
    female = fmt_count_pct(((baseline['sex'].astype(str).str.upper().isin(['F','FEMALE'])) | 
                        (baseline['sex'] == 1)).sum())
    edu = fmt_mean_sd(baseline['edu'])
    # apoe = fmt_count_pct((baseline['APOE4'].astype(str).str.upper().isin(['1','YES','TRUE',])).sum())
    apoe = fmt_count_pct(((baseline['APOE4'].astype(str).str.upper().isin(['1','YES','TRUE'])) | 
                        (baseline['APOE4'] == 1)).sum())
    # visits = fmt_iqr(df.groupby('id').size())
    visits = fmt_minmax(df.groupby('id').size())
    # followup = fmt_minmax(df.groupby('id')['months_since_baseline'].max() - df.groupby('id')['months_since_baseline'].min())
    followup = fmt_mean_sd(df.groupby('id')['months_since_baseline'].agg(lambda s: s.max(skipna=True) - s.min(skipna=True)).dropna())

    cu = fmt_count_pct((baseline['status'].astype(str).str.upper() == 'HC').sum())
    mci = fmt_count_pct((baseline['status'].astype(str).str.upper() == 'MCI').sum())
    ad = fmt_count_pct((baseline['status'].astype(str).str.upper() == 'AD').sum())
    meds = fmt_mean_sd(baseline['Total_Meds'])

    return pd.DataFrame({
        "Measure": [
            "Age at baseline (year)",
            "Gender (Female)",
            "Education (year)",
            "APOE4 (YES)**",
            "Visits (Record)",
            "Follow up intervals (Month)",
            "CU at baseline",
            "MCI at baseline",
            "AD at baseline",
            "Average medication taken at baseline"
        ],
        "Value": [
            age, female, edu, apoe, visits, followup, cu, mci, ad, meds
        ]
    })

# def baseline_age_by_med(df, med_cols, age_col='age', id_col='id', visit_col='visit_no'):
#     # 1) Keep baseline rows and ensure one row per participant
#     base = (df.loc[df[visit_col].eq(1), [id_col, age_col] + med_cols]
#               .sort_values(id_col)
#               .drop_duplicates(subset=id_col, keep='first'))
    
#     # 2) Reshape meds to long format: one row per id × drug
#     long = base.melt(id_vars=[id_col, age_col],
#                      value_name='user', var_name='drug')
    
#     # 3) Clean types: user → 0/1, drop missing ages or user flags
#     long['user'] = pd.to_numeric(long['user'], errors='coerce').astype('Int64')
#     long = long.dropna(subset=[age_col, 'user'])
    
#     # 4) Group and aggregate
#     out = (long.groupby(['drug', 'user'])[age_col]
#                .agg(n='count', mean='mean', sd='std')
#                .reset_index())
    
#     # 5) Tidy labels + formatting
#     out['user'] = out['user'].map({0: 'non-user', 1: 'user'})
#     out[['mean', 'sd']] = out[['mean', 'sd']].round(2)
#     out['mean_sd'] = out.apply(lambda r: f"{r['mean']} ({r['sd']})", axis=1)
    
#     # Optional: order columns
#     return out[['drug', 'user', 'n', 'mean', 'sd', 'mean_sd']]

import numpy as np
from scipy.stats import ttest_ind

# def baseline_age_by_med(df, med_cols, age_col='age', id_col='id', visit_col='visit_no'):
#     # 1) Keep baseline rows and ensure one row per participant
#     base = (df.loc[df[visit_col].eq(1), [id_col, age_col] + med_cols]
#               .sort_values(id_col)
#               .drop_duplicates(subset=id_col, keep='first'))

#     # 2) Reshape meds to long format: one row per id × drug
#     long = base.melt(id_vars=[id_col, age_col], value_name='user', var_name='drug')

#     # 3) Clean types: user → 0/1, drop missing ages or user flags
#     long['user'] = pd.to_numeric(long['user'], errors='coerce').astype('Int64')
#     long = long.dropna(subset=[age_col, 'user'])

#     # 4) Group and aggregate
#     out = (long.groupby(['drug', 'user'])[age_col]
#                .agg(n='count', mean='mean', sd='std')
#                .reset_index())

#     # 5) Tidy labels + formatting
#     out['user'] = out['user'].map({0: 'non-user', 1: 'user'})
#     out[['mean', 'sd']] = out[['mean', 'sd']].round(2)
#     out['mean_sd'] = out.apply(lambda r: f"{r['mean']} ({r['sd']})", axis=1)

#     # 6) Per-drug Welch t-test p-values (users vs non-users)
#     def welch_p(g):
#         a = g.loc[g['user'].eq(1), age_col].astype(float).dropna()
#         b = g.loc[g['user'].eq(0), age_col].astype(float).dropna()
#         if len(a) == 0 or len(b) == 0:
#             return np.nan
#         return ttest_ind(a, b, equal_var=False, nan_policy='omit').pvalue

#     pvals = (long.groupby('drug', as_index=False)
#                  .apply(lambda g: pd.Series({'p': welch_p(g)}))
#                  .reset_index(drop=True))

#     # 7) Merge and return ordered columns
#     out = out.merge(pvals, on='drug', how='left')
#     return out[['drug', 'user', 'n', 'mean', 'sd', 'mean_sd', 'p']]

def baseline_age_by_med(df, med_cols, age_col='age', id_col='id', visit_col='visit_no'):
    base = (df.loc[df[visit_col].eq(1), [id_col, age_col] + med_cols]
              .sort_values(id_col).drop_duplicates(subset=id_col, keep='first'))
    long = base.melt(id_vars=[id_col, age_col], value_name='user', var_name='drug')
    long['user'] = pd.to_numeric(long['user'], errors='coerce').astype('Int64')
    long = long.dropna(subset=[age_col, 'user'])

    out = (long.groupby(['drug', 'user'])[age_col]
               .agg(n='count', mean='mean', sd='std').reset_index())
    out['user'] = out['user'].map({0:'non-user',1:'user'})
    out[['mean','sd']] = out[['mean','sd']].round(2)
    out['mean_sd'] = out.apply(lambda r: f"{r['mean']} ({r['sd']})", axis=1)

    # Wald test per drug: H0 beta_user = 0 in OLS(age ~ 1 + user)
    def wald(g):
        g = g.dropna(subset=[age_col,'user']).copy()
        if g['user'].nunique() < 2: return pd.Series({'wald_chi2':np.nan,'wald_p':np.nan})
        X = sm.add_constant(g['user'].astype(float))
        y = g[age_col].astype(float)
        fit = sm.OLS(y, X, missing='drop').fit()
        w = fit.wald_test('user = 0')
        return pd.Series({'wald_chi2': float(w.statistic), 'wald_p': float(w.pvalue)})

    wald_df = long.groupby('drug', as_index=False).apply(wald).reset_index(drop=True)
    out = out.merge(wald_df, on='drug', how='left')
    return out[['drug','user','n','mean','sd','mean_sd','wald_chi2','wald_p']]


def qc(df, col, id_col='id', visit_col='visit_no', baseline_no=1):
    total_rows = len(df)
    missing = df[col].isna().sum()
    total_ids = df[id_col].nunique()

    has_1 = df.groupby(id_col)[col].apply(lambda s: s.notna().any()).sum()
    has_2 = df.groupby(id_col)[col].apply(lambda s: s.notna().sum() >= 2).sum()

    base_mask = df[visit_col].eq(baseline_no)
    base_ids_total = df.loc[base_mask, id_col].nunique()
    base_ids_with = df.loc[base_mask & df[col].notna(), id_col].nunique()

    print(f"Missing in {col}: {missing} / {total_rows}")
    print(f"Participants with ≥1 '{col}' value: {has_1} / {total_ids}")
    print(f"Participants with ≥2 '{col}' values: {has_2} / {total_ids}")
    print(f"Participants with baseline '{col}' value: {base_ids_with} / {base_ids_total}")