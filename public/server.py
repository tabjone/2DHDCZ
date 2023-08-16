from flask import Flask, request
import subprocess
import os

app = Flask(__name__)

@app.route('/run', methods=['POST'])
def run():
    # Get the filename and time from the request
    filename = request.form.get('filename')
    time = request.form.get('time')

    # Ensure the filename and time are provided
    if filename is None or time is None:
        return 'Error: Missing filename or time.', 400

    # Ensure the file exists
    if not os.path.exists(filename):
        return 'Error: File does not exist.', 400

    # Ensure the time input is a float
    try:
        float(time)
    except ValueError:
        return 'Error: Invalid time value.', 400

    # Run the command
    result = subprocess.run(['./main.o', filename, time], capture_output=True, text=True)

    if result.returncode != 0:
        return f'Error: {result.stderr}', 500

    return result.stdout

if __name__ == '__main__':
    app.run(port=8080)
