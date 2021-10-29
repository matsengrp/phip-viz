FROM quay.io/matsengrp/python3.7
WORKDIR /app
RUN python -m pip install --upgrade pip
COPY requirements.txt ./requirements.txt
RUN python -m pip install -r requirements.txt
RUN python -m pip install git+https://github.com/matsengrp/phippery.git@416406d693aaf0c5d15893fb2755e3c66648c8c0 
EXPOSE 8501
COPY streamlit_app.py ./streamlit_app.py
ENTRYPOINT ["streamlit", "run"]
CMD ["streamlit_app.py"]
