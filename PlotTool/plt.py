"""
    plot for: 
    ---
        - bar
        - pie
"""

from __future__ import annotations
from pandas.core.frame import DataFrame
import pandas as pd
from   matplotlib import pyplot as plt
import seaborn as sns
from   typing import List, Dict, Tuple, Union, Optional, Any

def plot_bar(data: DataFrame, title: str, x: str, y: str, x_label: str, y_label: str, percentage: bool = False) -> None:
    """
        plot_bar:
        ---
            input:
        ---
                data : data (DataFrame)
                title : plot_title_name (str)
                x : DataFrame column name for x axis (str)
                y : DataFrame column name for y axis (str)
                x_label : x_label_name (str)
                y_label : y_label_name (str)
                percentage : Whether to display the percentage of each bar relative to the total (bool)
        ---
            output: plot
        ---
    """
    plt.figure(figsize=(10, 6))
    # sns.barplot(x=x, y=y, data=data, palette="viridis")
    ax = sns.barplot(x=x, y=y, data=data, palette="viridis")
    plt.xticks(rotation=45)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    # NOTE: 添加百分比 (根據參數是否設定 precentage 為 True 決定)
    # If percentage is True, annotate each bar with its percentage of the total
    if percentage:
        total = sum(data[y])
        for p in ax.patches:
            height = p.get_height()
            ax.text(p.get_x() + p.get_width() / 2.,
                    height + 3,
                    '{:1.2f}%'.format(100 * height / total),
                    ha="center")

    plt.tight_layout()
    plt.show()

def plot_pie(data: DataFrame, title: str, x: str, y: str, x_label: str, y_label: str) -> None:
    """
        plot_pie:
        ---
            input:
        ---
                data : data (DataFrame)
                title : plot_title_name (str)
                x : DataFrame column name for pie labels (str)
                y : DataFrame column name for pie values (str)
                x_label : Legend title (str)
                y_label : Not used in pie chart but kept for consistency (str)
        ---
            output: plot
        ---
    """
    plt.figure(figsize=(10, 6))
    colors = sns.color_palette("viridis", len(data[y]))
    wedges, texts, autotexts = plt.pie(data[y], labels=data[x], autopct='%1.1f%%', startangle=140, colors=colors)
    
    # Add legend
    plt.legend(wedges, data[x], title=x_label, loc="upper right", bbox_to_anchor=(1, 0, 0.5, 1))
    
    plt.title(title)
    plt.tight_layout()
    plt.show()