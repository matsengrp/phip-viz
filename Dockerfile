FROM quay.io/matsengrp/python3.7
WORKDIR /app
RUN python -m pip install --upgrade pip
COPY requirements.txt ./requirements.txt
RUN python -m pip install git+https://github.com/matsengrp/phippery.git
RUN python -m pip install -r requirements.txt
EXPOSE 8501
COPY streamlit_app.py ./streamlit_app.py
ENTRYPOINT ["streamlit", "run"]
CMD ["streamlit_app.py"]
