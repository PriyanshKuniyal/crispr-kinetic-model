import os
import sys

def try_extract():
    # Try PyMuPDF (fitz) first, then pdfplumber, then pypdf, then PyPDF2
    methods = []
    
    # Method 1: fitz
    try:
        import fitz
        methods.append(('fitz', fitz))
        print("fitz available")
    except ImportError:
        pass
        
    # Method 2: pdfplumber
    try:
        import pdfplumber
        methods.append(('pdfplumber', pdfplumber))
        print("pdfplumber available")
    except ImportError:
        pass
        
    # Method 3: pypdf
    try:
        import pypdf
        methods.append(('pypdf', pypdf))
        print("pypdf available")
    except ImportError:
        pass

    # Method 4: PyPDF2
    try:
        import PyPDF2
        methods.append(('PyPDF2', PyPDF2))
        print("PyPDF2 available")
    except ImportError:
        pass

    if not methods:
        print("No PDF reading libraries found. Trying to install pypdf...")
        import subprocess
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pypdf"])
        import pypdf
        methods.append(('pypdf', pypdf))

    name, lib = methods[0]
    print(f"Using {name} to extract text.")
    
    pdf_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "nature_communications_paper.pdf")
    
    text_content = []
    if name == 'fitz':
        doc = lib.open(pdf_path)
        for i, page in enumerate(doc):
            text_content.append(f"--- PAGE {i+1} ---")
            text_content.append(page.get_text())
    elif name == 'pdfplumber':
        with lib.open(pdf_path) as pdf:
            for i, page in enumerate(pdf.pages):
                text_content.append(f"--- PAGE {i+1} ---")
                text_content.append(page.extract_text())
    elif name == 'pypdf':
        reader = lib.PdfReader(pdf_path)
        for i, page in enumerate(reader.pages):
            text_content.append(f"--- PAGE {i+1} ---")
            text_content.append(page.extract_text())
    elif name == 'PyPDF2':
        reader = lib.PdfReader(pdf_path)
        for i, page in enumerate(reader.pages):
            text_content.append(f"--- PAGE {i+1} ---")
            text_content.append(page.extract_text())

    full_text = "\n".join(text_content)
    with open("paper_text.txt", "w", encoding="utf-8") as f:
        f.write(full_text)
    print(f"Extracted {len(full_text)} characters and wrote to paper_text.txt")

if __name__ == "__main__":
    try_extract()
