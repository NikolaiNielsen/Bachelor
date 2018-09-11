import sys
from PyQt5 import QtWidgets, QtCore


class window(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setGeometry(50, 50, 500, 300)
        self.setWindowTitle("yoyoy")
        self.home()

    def home(self):
        btn = QtWidgets.QPushButton("Quitter", self)
        btn.clicked.connect(QtCore.QCoreApplication.instance().quit)
        btn.resize(100, 100)
        btn.move(100, 100)
        self.show()


def main():
    app = QtWidgets.QApplication(sys.argv)
    GUI = window()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
