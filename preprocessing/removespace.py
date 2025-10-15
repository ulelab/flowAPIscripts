#!/usr/bin/env python3

import sys
import subprocess

def process_fastq(input_file, output_file):
    """Process FASTQ file to remove spaces from headers using pigz."""
    is_gzipped = input_file.endswith('.gz')
    
    # Set up input stream
    if is_gzipped:
        pigz_decompress = subprocess.Popen(['pigz', '-dc', input_file], 
                                         stdout=subprocess.PIPE)
    else:
        in_file = open(input_file, 'rb')
        pigz_decompress = subprocess.Popen(['cat'], 
                                         stdin=in_file,
                                         stdout=subprocess.PIPE)

    # Set up output stream
    out_file_name = output_file if output_file else input_file.replace('.fastq', '.cleaned.fastq')
    if not out_file_name.endswith('.gz'):
        out_file_name += '.gz'
        
    pigz_compress = subprocess.Popen(['pigz', '-c'], 
                                   stdin=subprocess.PIPE, 
                                   stdout=open(out_file_name, 'wb'))
    
    try:
        line_count = 0
        for line in pigz_decompress.stdout:
            if line_count % 4 == 0:  # Header line
                cleaned_header = line.decode().strip().replace(' ', '_').encode() + b'\n'
                pigz_compress.stdin.write(cleaned_header)
            else:
                pigz_compress.stdin.write(line)
            line_count += 1
            
    finally:
        pigz_decompress.stdout.close()
        pigz_decompress.wait()
        if not is_gzipped:
            in_file.close()
        
        pigz_compress.stdin.close()
        pigz_compress.wait()

def main():
    if len(sys.argv) < 2:
        print("Usage: python script.py input.fastq [output.fastq]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    try:
        process_fastq(input_file, output_file)
    except Exception as e:
        print(f"Error processing file: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()