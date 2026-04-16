#!/usr/bin/env python3
"""
Extract tables from Randal's paper to compare with our results.
"""

import PyPDF2
import re

pdf_path = '/home/por07g/Documents/Projects/GBR_modeling/Connectivity/R_codes/Randal/s42003-024-05824-3.pdf'

# Extract text from PDF
with open(pdf_path, 'rb') as file:
    pdf_reader = PyPDF2.PdfReader(file)
    
    print(f"Total pages: {len(pdf_reader.pages)}\n")
    
    # Search for tables with LD50 or competency data
    for page_num in range(len(pdf_reader.pages)):
        page = pdf_reader.pages[page_num]
        text = page.extract_text()
        
        # Look for keywords that might indicate the table
        if any(word in text.lower() for word in ['ld50', 'competenc', 'settlement', 'median', 'table']):
            print(f"\n{'='*70}")
            print(f"PAGE {page_num + 1}")
            print(f"{'='*70}")
            print(text)
            print("\n")

